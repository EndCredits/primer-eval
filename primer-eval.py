#!/usr/bin/env python3
import re
from dataclasses import dataclass, asdict
from typing import List, Tuple, Optional
from datetime import datetime
from primer3 import (
    calc_tm, calc_hairpin, calc_homodimer,
    calc_heterodimer, calc_end_stability
)
from argparse import ArgumentParser

# Convert cal to kcal
# Primer3-py returns cal by default
CAL_TO_KCAL = 1000.0  # 1 kcal = 1000 cal


@dataclass
class PrimerAnalysis:
    """Primer analysis result"""
    sequence: str
    length: int
    gc_content: float
    tm: float
    hairpin_dg: float  # kcal/mol
    homodimer_dg: float  # kcal/mol

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class PrimerPairAnalysis:
    """引物对的分析结果"""
    forward: PrimerAnalysis
    reverse: PrimerAnalysis
    tm_difference: float
    heterodimer_dg: float  # kcal/mol
    three_prime_risk: bool
    three_prime_details: List[str]
    warnings: List[str]
    passed: bool

    def to_dict(self) -> dict:
        return {
            "forward": self.forward.to_dict(),
            "reverse": self.reverse.to_dict(),
            "pair": {
                "tm_difference": self.tm_difference,
                "heterodimer_dg": self.heterodimer_dg,
                "three_prime_risk": self.three_prime_risk,
                "three_prime_details": self.three_prime_details
            },
            "warnings": self.warnings,
            "passed": self.passed
        }


class Primer3Validator:
    """
    A simple primer evaluation tool based on Primer3-Py.
    """

    def __init__(self, **kwargs):
        # 默认参数（符合Primer3标准）
        self.params = {
            'mv_conc': 50.0,    # Monovalent ion (mM)
            'dv_conc': 1.5,     # Divalent ion (mM)
            'dntp_conc': 0.6,   # dNTP concentration (mM)
            'dna_conc': 50.0,   # Primer concentration (nM)
            'temp_c': 37.0,     # Temperature for ΔG calculations (°C)
            'max_loop': 30,     # Maximum loop size for secondary structures
            'output_structure': True,  # Need structure info for ΔG
            'max_nn_length': 60,       # Max length for nearest-neighbor calcs
            'tm_method': 'santalucia',  # Tm calculation method
            'salt_corrections_method': 'santalucia'  # Salt correction method
        }

        # QC thresholds (kcal/mol)
        self.thresholds = {
            'min_gc': 40.0,      # Minimum GC%
            'max_gc': 60.0,      # Maximum GC%
            'max_tm_diff': 5.0,  # Maximum Tm difference between primers
            'max_hairpin_dg': -3.0,    # Hairpin stability threshold
            'max_homodimer_dg': -5.0,  # Homodimer stability threshold
            'max_heterodimer_dg': -5.0,  # Heterodimer stability threshold
            'max_three_prime_dg': -3.0  # 3'-end stability threshold
        }

        # Allows overriding default parameters and thresholds
        for key, value in kwargs.items():
            if key in self.params:
                self.params[key] = value
            elif key in self.thresholds:
                self.thresholds[key] = value

    def _validate_sequence(self, seq: str, name: str = "Primer") -> str:
        """Validate and normalize DNA sequences"""
        if not seq:
            raise ValueError(f"{name} sequence cannot be empty")

        seq = seq.upper().strip()

        # Check effective bases
        valid_bases = set("ATGCNRYKMSWBDHV")
        if not all(base in valid_bases for base in seq):
            invalid_bases = set(seq) - valid_bases
            raise ValueError(f"{name} contains invalid bases: {invalid_bases}")

        # Check the length limit, as Primer3 requires <60bp for thermodynamic calculations.
        if len(seq) >= 60:
            raise ValueError(
                f"{name} sequence must be < 60 bp for reliable thermodynamic calculations (current: {len(seq)} bp)")
        if len(seq) < 15:
            raise ValueError(
                f"{name} sequence should be at least 15 bp for PCR applications (current: {len(seq)} bp)")

        return seq

    def _calculate_gc_content(self, seq: str) -> float:
        """Calculate GC content and handle degenerate bases."""
        gc_count = 0
        total = len(seq)

        for base in seq:
            if base in 'GC':  # G or C
                gc_count += 1
            elif base == 'S':  # S = G or C
                gc_count += 1
            elif base == 'N':  # N = A, T, G, or C (50% chance of G/C)
                gc_count += 0.5

        return (gc_count / total) * 100 if total > 0 else 0

    def _convert_cal_to_kcal(self, value: float) -> float:
        return value / CAL_TO_KCAL

    def _safe_thermo_calculation(self, func, *args, **kwargs) -> Tuple[float, bool]:
        try:
            result = func(*args, **kwargs)
            if hasattr(result, 'structure_found') and result.structure_found:
                if hasattr(result, 'dg'):
                    return self._convert_cal_to_kcal(result.dg), True
            return 0.0, False
        except RuntimeError as e:
            if "sequence too long" in str(e).lower():
                raise ValueError(
                    f"Sequence too long for thermodynamic analysis: {args[0]}")
            return 0.0, False
        except Exception as e:
            return 0.0, False

    def _analyze_single_primer(self, seq: str) -> PrimerAnalysis:
        # Tm calculation (using different parameters just for Tm calculation)
        tm_params = {
            'mv_conc': self.params['mv_conc'],
            'dv_conc': self.params['dv_conc'],
            'dntp_conc': self.params['dntp_conc'],
            'dna_conc': self.params['dna_conc'],
            'max_nn_length': self.params['max_nn_length'],
            'tm_method': self.params['tm_method'],
            'salt_corrections_method': self.params['salt_corrections_method']
        }

        try:
            tm = calc_tm(seq, **tm_params)
        except Exception as e:
            raise RuntimeError(
                f"Tm calculation failed for sequence '{seq}': {str(e)}")

        # Thermodynamic calculation parameters just for hairpins, dimers, etc.
        thermo_params = {
            'mv_conc': self.params['mv_conc'],
            'dv_conc': self.params['dv_conc'],
            'dntp_conc': self.params['dntp_conc'],
            'dna_conc': self.params['dna_conc'],
            'temp_c': self.params['temp_c'],
            'max_loop': self.params['max_loop'],
            'output_structure': self.params['output_structure']
        }

        # harpin
        hairpin_dg, hairpin_found = self._safe_thermo_calculation(
            calc_hairpin, seq, **thermo_params)

        # homodimers
        homodimer_dg, homodimer_found = self._safe_thermo_calculation(
            calc_homodimer, seq, **thermo_params)

        return PrimerAnalysis(
            sequence=seq,
            length=len(seq),
            gc_content=round(self._calculate_gc_content(seq), 2),
            tm=round(tm, 2),
            hairpin_dg=round(hairpin_dg, 2),
            homodimer_dg=round(homodimer_dg, 2)
        )

    def _check_three_prime_stability(self, seq1: str, seq2: str) -> Tuple[bool, float, str]:
        """Check the stability of the 3' end"""
        thermo_params = {
            'mv_conc': self.params['mv_conc'],
            'dv_conc': self.params['dv_conc'],
            'dntp_conc': self.params['dntp_conc'],
            'dna_conc': self.params['dna_conc'],
            'temp_c': self.params['temp_c'],
            'max_loop': self.params['max_loop']
        }

        try:
            result = calc_end_stability(seq1, seq2, **thermo_params)
            if hasattr(result, 'structure_found') and result.structure_found and hasattr(result, 'dg'):
                dg_kcal = self._convert_cal_to_kcal(result.dg)
                is_risky = dg_kcal < self.thresholds['max_three_prime_dg']
                detail = f"3'-end stability ΔG = {dg_kcal:.2f} kcal/mol"
                return is_risky, dg_kcal, detail
        except Exception as e:
            pass

        return False, 0.0, ""

    def analyze_primer_pair(self, forward_seq: str, reverse_seq: str) -> PrimerPairAnalysis:
        fwd_seq = self._validate_sequence(forward_seq, "Forward primer")
        rev_seq = self._validate_sequence(reverse_seq, "Reverse primer")

        warnings = []

        fwd_analysis = self._analyze_single_primer(fwd_seq)
        rev_analysis = self._analyze_single_primer(rev_seq)

        thermo_params = {
            'mv_conc': self.params['mv_conc'],
            'dv_conc': self.params['dv_conc'],
            'dntp_conc': self.params['dntp_conc'],
            'dna_conc': self.params['dna_conc'],
            'temp_c': self.params['temp_c'],
            'max_loop': self.params['max_loop'],
            'output_structure': self.params['output_structure']
        }

        heterodimer_dg, heterodimer_found = self._safe_thermo_calculation(
            calc_heterodimer, fwd_seq, rev_seq, **thermo_params
        )

        fwd_to_rev_risk, fwd_dg, fwd_detail = self._check_three_prime_stability(
            fwd_seq, rev_seq)
        rev_to_fwd_risk, rev_dg, rev_detail = self._check_three_prime_stability(
            rev_seq, fwd_seq)

        three_prime_risk = fwd_to_rev_risk or rev_to_fwd_risk
        three_prime_details = []
        if fwd_to_rev_risk:
            three_prime_details.append(f"Forward 3'-end: {fwd_detail}")
        if rev_to_fwd_risk:
            three_prime_details.append(f"Reverse 3'-end: {rev_detail}")

        tm_diff = abs(fwd_analysis.tm - rev_analysis.tm)

        if not (self.thresholds['min_gc'] <= fwd_analysis.gc_content <= self.thresholds['max_gc']):
            warnings.append(f"[GC%] Forward primer: {fwd_analysis.gc_content:.1f}% "
                            f"(recommended: {self.thresholds['min_gc']}-{self.thresholds['max_gc']}%)")

        if not (self.thresholds['min_gc'] <= rev_analysis.gc_content <= self.thresholds['max_gc']):
            warnings.append(f"[GC%] Reverse primer: {rev_analysis.gc_content:.1f}% "
                            f"(recommended: {self.thresholds['min_gc']}-{self.thresholds['max_gc']}%)")

        if tm_diff > self.thresholds['max_tm_diff']:
            warnings.append(f"[Tm] Tm difference: {tm_diff:.1f}°C "
                            f"(maximum allowed: {self.thresholds['max_tm_diff']}°C)")

        if fwd_analysis.hairpin_dg < self.thresholds['max_hairpin_dg']:
            warnings.append(f"[Hairpin] Forward primer: ΔG = {fwd_analysis.hairpin_dg:.2f} kcal/mol "
                            f"(threshold: {self.thresholds['max_hairpin_dg']} kcal/mol)")

        if rev_analysis.hairpin_dg < self.thresholds['max_hairpin_dg']:
            warnings.append(f"[Hairpin] Reverse primer: ΔG = {rev_analysis.hairpin_dg:.2f} kcal/mol "
                            f"(threshold: {self.thresholds['max_hairpin_dg']} kcal/mol)")

        if fwd_analysis.homodimer_dg < self.thresholds['max_homodimer_dg']:
            warnings.append(f"[Homodimer] Forward primer: ΔG = {fwd_analysis.homodimer_dg:.2f} kcal/mol "
                            f"(threshold: {self.thresholds['max_homodimer_dg']} kcal/mol)")

        if rev_analysis.homodimer_dg < self.thresholds['max_homodimer_dg']:
            warnings.append(f"[Homodimer] Reverse primer: ΔG = {rev_analysis.homodimer_dg:.2f} kcal/mol "
                            f"(threshold: {self.thresholds['max_homodimer_dg']} kcal/mol)")

        if heterodimer_dg < self.thresholds['max_heterodimer_dg']:
            warnings.append(f"[Heterodimer] Primer pair: ΔG = {heterodimer_dg:.2f} kcal/mol "
                            f"(threshold: {self.thresholds['max_heterodimer_dg']} kcal/mol)")

        if three_prime_risk:
            for detail in three_prime_details:
                warnings.append(f"[3'-end risk] {detail}")

        passed = len(warnings) == 0

        return PrimerPairAnalysis(
            forward=fwd_analysis,
            reverse=rev_analysis,
            tm_difference=round(tm_diff, 2),
            heterodimer_dg=round(heterodimer_dg, 2),
            three_prime_risk=three_prime_risk,
            three_prime_details=three_prime_details,
            warnings=warnings,
            passed=passed
        )

    def generate_report(self, analysis: PrimerPairAnalysis) -> str:
        report_lines = []

        report_lines.append("=" * 60)
        report_lines.append("PRIMER PAIR QUALITY ASSESSMENT REPORT")
        report_lines.append("=" * 60)
        report_lines.append(
            f"Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append("-" * 60)

        report_lines.append("\nFORWARD PRIMER")
        report_lines.append("-" * 30)
        report_lines.append(f"Sequence: 5'-{analysis.forward.sequence}-3'")
        report_lines.append(f"Length:   {analysis.forward.length} bp")
        report_lines.append(f"GC%:      {analysis.forward.gc_content:.2f}%")
        report_lines.append(f"Tm:       {analysis.forward.tm:.2f} °C")
        report_lines.append(
            f"Hairpin ΔG:  {analysis.forward.hairpin_dg:.2f} kcal/mol")
        report_lines.append(
            f"Homodimer ΔG: {analysis.forward.homodimer_dg:.2f} kcal/mol")

        report_lines.append("\nREVERSE PRIMER")
        report_lines.append("-" * 30)
        report_lines.append(f"Sequence: 5'-{analysis.reverse.sequence}-3'")
        report_lines.append(f"Length:   {analysis.reverse.length} bp")
        report_lines.append(f"GC%:      {analysis.reverse.gc_content:.2f}%")
        report_lines.append(f"Tm:       {analysis.reverse.tm:.2f} °C")
        report_lines.append(
            f"Hairpin ΔG:  {analysis.reverse.hairpin_dg:.2f} kcal/mol")
        report_lines.append(
            f"Homodimer ΔG: {analysis.reverse.homodimer_dg:.2f} kcal/mol")

        report_lines.append("\nPRIMER PAIR")
        report_lines.append("-" * 30)
        report_lines.append(
            f"Tm difference:    {analysis.tm_difference:.2f} °C")
        report_lines.append(
            f"Heterodimer ΔG:   {analysis.heterodimer_dg:.2f} kcal/mol")

        if analysis.three_prime_risk:
            report_lines.append("3'-end stability risks detected:")
            for detail in analysis.three_prime_details:
                report_lines.append(f"  • {detail}")
        else:
            report_lines.append(
                "3'-end stability: No significant risks detected")

        report_lines.append("\nASSESSMENT RESULT")
        report_lines.append("-" * 30)

        if analysis.passed:
            report_lines.append("✅ PASSED - No significant issues detected")
        else:
            report_lines.append("❌ FAILED - Issues detected:")
            for i, warning in enumerate(analysis.warnings, 1):
                report_lines.append(f"  {i}. {warning}")

        report_lines.append("=" * 60)
        return "\n".join(report_lines)


# Testcase
if __name__ == "__main__":
    validator = Primer3Validator(
        mv_conc=50.0,    # Monovalent ions (Na+/K+)
        dv_conc=1.5,     # Divalent ions (Mg2+)
        dntp_conc=0.6,   # dNTP concentration
        dna_conc=50.0    # Primer concentration (nM)
    )

    argparser = ArgumentParser(
        "primer-eval", description="A simple tool to evaluation primer properites")
    argparser.add_argument(
        "forward_primer", help="Sequence of your forward primer (5\'-3\')")
    argparser.add_argument(
        "reverse_primer", help="Sequence of your reverse primer (5\'-3\')")
    args = argparser.parse_args()

    forward_primer = str(args.forward_primer).upper()
    reverse_primer = str(args.reverse_primer).upper()

    try:
        analysis = validator.analyze_primer_pair(
            forward_primer, reverse_primer)

        report = validator.generate_report(analysis)
        print(report)

        json_result = analysis.to_dict()
        print("\n✅ Analysis completed successfully!")
        print(
            f"✅ JSON result available with {len(json_result['warnings'])} warnings")

    except ValueError as e:
        print(f"❌ Validation error: {e}")
    except RuntimeError as e:
        print(f"❌ Analysis error: {e}")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
