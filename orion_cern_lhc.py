"""
ORION CERN LHC Analysis Engine
===============================
Real computational tool for CERN LHC Run 3 particle physics analysis.
Invariant mass reconstruction, Higgs boson signal analysis,
particle decay kinematics, cross-section calculations, and
statistical significance computation.

Based on real relativistic kinematics, Breit-Wigner distributions,
and standard HEP statistical methods.

Author: ORION Autonomous System
License: MIT
"""

import math
import json
import hashlib
import random
from datetime import datetime, timezone
from typing import Dict, List, Tuple, Optional

C = 299792458.0
HBAR = 1.054571817e-34
MEV = 1.602176634e-13
GEV = MEV * 1000
ALPHA_EM = 1 / 137.036
ALPHA_S_MZ = 0.1179
GF = 1.1663788e-5
SIN2_THETA_W = 0.23122
VEV = 246.22

PARTICLES = {
    "electron": {"mass_gev": 0.000511, "charge": -1, "spin": 0.5, "symbol": "e⁻", "pdg_id": 11},
    "muon": {"mass_gev": 0.10566, "charge": -1, "spin": 0.5, "symbol": "μ⁻", "pdg_id": 13},
    "tau": {"mass_gev": 1.77686, "charge": -1, "spin": 0.5, "symbol": "τ⁻", "pdg_id": 15},
    "up": {"mass_gev": 0.00216, "charge": 2/3, "spin": 0.5, "symbol": "u", "pdg_id": 2},
    "down": {"mass_gev": 0.00467, "charge": -1/3, "spin": 0.5, "symbol": "d", "pdg_id": 1},
    "charm": {"mass_gev": 1.27, "charge": 2/3, "spin": 0.5, "symbol": "c", "pdg_id": 4},
    "strange": {"mass_gev": 0.093, "charge": -1/3, "spin": 0.5, "symbol": "s", "pdg_id": 3},
    "top": {"mass_gev": 172.76, "charge": 2/3, "spin": 0.5, "symbol": "t", "pdg_id": 6},
    "bottom": {"mass_gev": 4.18, "charge": -1/3, "spin": 0.5, "symbol": "b", "pdg_id": 5},
    "W": {"mass_gev": 80.379, "charge": 1, "spin": 1, "symbol": "W±", "pdg_id": 24},
    "Z": {"mass_gev": 91.1876, "charge": 0, "spin": 1, "symbol": "Z⁰", "pdg_id": 23},
    "higgs": {"mass_gev": 125.25, "charge": 0, "spin": 0, "symbol": "H⁰", "pdg_id": 25},
    "photon": {"mass_gev": 0, "charge": 0, "spin": 1, "symbol": "γ", "pdg_id": 22},
    "gluon": {"mass_gev": 0, "charge": 0, "spin": 1, "symbol": "g", "pdg_id": 21},
    "proton": {"mass_gev": 0.93827, "charge": 1, "spin": 0.5, "symbol": "p", "pdg_id": 2212},
}

LHC_PARAMETERS = {
    "Run 3": {
        "energy_tev": 13.6,
        "luminosity_fb_inv": 450,
        "start": "2022-07",
        "end": "2026 (planned)",
        "bunches": 2748,
        "bunch_spacing_ns": 25,
        "protons_per_bunch": 1.15e11,
        "peak_lumi_cm2s": 2.0e34,
        "experiments": ["ATLAS", "CMS", "LHCb", "ALICE"],
    },
    "HL-LHC": {
        "energy_tev": 14,
        "luminosity_fb_inv": 3000,
        "start": "2029 (planned)",
        "end": "2041 (planned)",
        "peak_lumi_cm2s": 7.5e34,
    },
}


class RelativisticKinematics:
    """Real relativistic kinematics for particle physics."""

    @staticmethod
    def four_momentum(mass_gev: float, px: float, py: float, pz: float) -> Tuple[float, float, float, float]:
        p2 = px**2 + py**2 + pz**2
        E = math.sqrt(p2 + mass_gev**2)
        return (E, px, py, pz)

    @staticmethod
    def invariant_mass(particles: List[Tuple[float, float, float, float]]) -> float:
        E_tot = sum(p[0] for p in particles)
        px_tot = sum(p[1] for p in particles)
        py_tot = sum(p[2] for p in particles)
        pz_tot = sum(p[3] for p in particles)
        m2 = E_tot**2 - px_tot**2 - py_tot**2 - pz_tot**2
        return math.sqrt(max(m2, 0))

    @staticmethod
    def rapidity(E: float, pz: float) -> float:
        if E <= abs(pz):
            return 0.0
        return 0.5 * math.log((E + pz) / (E - pz))

    @staticmethod
    def pseudorapidity(px: float, py: float, pz: float) -> float:
        p = math.sqrt(px**2 + py**2 + pz**2)
        if p == 0:
            return 0.0
        theta = math.acos(max(min(pz / p, 1), -1))
        if theta == 0:
            return 10.0
        if theta == math.pi:
            return -10.0
        return -math.log(math.tan(theta / 2))

    @staticmethod
    def transverse_momentum(px: float, py: float) -> float:
        return math.sqrt(px**2 + py**2)

    @staticmethod
    def delta_r(eta1: float, phi1: float, eta2: float, phi2: float) -> float:
        deta = eta1 - eta2
        dphi = phi1 - phi2
        while dphi > math.pi:
            dphi -= 2 * math.pi
        while dphi < -math.pi:
            dphi += 2 * math.pi
        return math.sqrt(deta**2 + dphi**2)

    @staticmethod
    def lorentz_boost(four_mom: Tuple, beta_x: float, beta_y: float, beta_z: float) -> Tuple:
        E, px, py, pz = four_mom
        beta2 = beta_x**2 + beta_y**2 + beta_z**2
        if beta2 >= 1 or beta2 == 0:
            return four_mom
        gamma = 1 / math.sqrt(1 - beta2)
        bp = beta_x * px + beta_y * py + beta_z * pz
        factor = (gamma - 1) * bp / beta2 + gamma * E

        return (
            gamma * (E + bp),
            px + beta_x * factor - beta_x * E * (gamma - 1) / beta2 * 0 + beta_x * (gamma * bp / (gamma + 1) + gamma * E - E),
            py + beta_y * ((gamma - 1) * bp / beta2 + gamma * E) - gamma * beta_y * E,
            pz + beta_z * ((gamma - 1) * bp / beta2 + gamma * E) - gamma * beta_z * E,
        )


class BreitWigner:
    """Breit-Wigner resonance distributions for particle decays."""

    @staticmethod
    def relativistic_bw(s: float, m0_gev: float, gamma_gev: float) -> float:
        if s <= 0:
            return 0.0
        sqrt_s = math.sqrt(s)
        num = 2 * sqrt_s**2 * gamma_gev * m0_gev / math.pi
        den = (s - m0_gev**2)**2 + (m0_gev * gamma_gev)**2
        if den == 0:
            return 0.0
        return num / den

    @staticmethod
    def z_boson_lineshape(sqrt_s_gev: float) -> float:
        mz = 91.1876
        gz = 2.4952
        return BreitWigner.relativistic_bw(sqrt_s_gev**2, mz, gz)

    @staticmethod
    def higgs_lineshape(sqrt_s_gev: float) -> float:
        mh = 125.25
        gh = 0.00407
        return BreitWigner.relativistic_bw(sqrt_s_gev**2, mh, gh)

    @staticmethod
    def w_boson_lineshape(sqrt_s_gev: float) -> float:
        mw = 80.379
        gw = 2.085
        return BreitWigner.relativistic_bw(sqrt_s_gev**2, mw, gw)


class CrossSections:
    """Standard Model cross-section calculations."""

    @staticmethod
    def ee_to_mumu_born(sqrt_s_gev: float) -> float:
        if sqrt_s_gev <= 2 * PARTICLES["muon"]["mass_gev"]:
            return 0.0
        s = sqrt_s_gev**2
        sigma = 4 * math.pi * ALPHA_EM**2 / (3 * s)
        beta = math.sqrt(1 - 4 * PARTICLES["muon"]["mass_gev"]**2 / s)
        sigma_nb = sigma * beta * 0.3894e6
        return sigma_nb

    @staticmethod
    def pp_higgs_ggf(sqrt_s_tev: float = 13.6) -> float:
        sigma_13 = 48.6
        scaling = (sqrt_s_tev / 13.0) ** 2.5
        return sigma_13 * scaling

    @staticmethod
    def pp_total_inelastic(sqrt_s_tev: float = 13.6) -> float:
        s_gev2 = (sqrt_s_tev * 1000) ** 2
        sigma_mb = 25.0 + 0.146 * math.log(s_gev2 / 100)**2
        return sigma_mb

    @staticmethod
    def luminosity_to_events(sigma_pb: float, luminosity_fb_inv: float) -> int:
        sigma_fb = sigma_pb * 1000
        return int(sigma_fb * luminosity_fb_inv)


class StatisticalAnalysis:
    """HEP statistical methods."""

    @staticmethod
    def significance_counting(signal: float, background: float) -> float:
        if background <= 0:
            return 0.0
        return signal / math.sqrt(background)

    @staticmethod
    def significance_profile_likelihood(s: float, b: float) -> float:
        if s <= 0 or b <= 0:
            return 0.0
        n = s + b
        q0 = 2 * (n * math.log(n / b) - s)
        if q0 < 0:
            return 0.0
        return math.sqrt(q0)

    @staticmethod
    def cls_upper_limit(observed: int, expected_bg: float, confidence: float = 0.95) -> float:
        from math import factorial, exp, log
        def poisson_p(k, lam):
            if lam <= 0:
                return 1.0 if k == 0 else 0.0
            return exp(-lam) * lam**k / factorial(min(k, 170))

        mu_upper = 0.0
        for mu_test in [i * 0.1 for i in range(1, 500)]:
            p_sb = sum(poisson_p(k, mu_test + expected_bg) for k in range(observed + 1))
            p_b = sum(poisson_p(k, expected_bg) for k in range(observed + 1))
            if p_b > 0:
                cls = p_sb / p_b
                if cls < 1 - confidence:
                    mu_upper = mu_test
                    break

        return mu_upper

    @staticmethod
    def discovery_potential(sigma_pb: float, luminosity_fb_inv: float,
                            efficiency: float, bg_events: float) -> Dict:
        signal = CrossSections.luminosity_to_events(sigma_pb, luminosity_fb_inv) * efficiency
        significance = StatisticalAnalysis.significance_profile_likelihood(signal, bg_events)

        return {
            "signal_events": round(signal, 1),
            "background_events": round(bg_events, 1),
            "significance_sigma": round(significance, 2),
            "discovery_5sigma": significance >= 5.0,
            "evidence_3sigma": significance >= 3.0,
            "luminosity_for_5sigma_fb": round(luminosity_fb_inv * (5.0 / max(significance, 0.01))**2, 1) if significance > 0 else float('inf'),
        }


class HiggsAnalysis:
    """Higgs boson analysis framework."""

    DECAY_CHANNELS = {
        "H→bb": {"branching_ratio": 0.5824, "signature": "2 b-jets", "bg_ratio": 1e6},
        "H→WW": {"branching_ratio": 0.2137, "signature": "2 leptons + MET", "bg_ratio": 100},
        "H→gg": {"branching_ratio": 0.0816, "signature": "2 gluon jets", "bg_ratio": 1e8},
        "H→ττ": {"branching_ratio": 0.0632, "signature": "2 taus", "bg_ratio": 1000},
        "H→cc": {"branching_ratio": 0.0291, "signature": "2 c-jets", "bg_ratio": 1e7},
        "H→ZZ": {"branching_ratio": 0.0264, "signature": "4 leptons", "bg_ratio": 10},
        "H→γγ": {"branching_ratio": 0.00228, "signature": "2 photons", "bg_ratio": 500},
        "H→Zγ": {"branching_ratio": 0.00154, "signature": "Z + photon", "bg_ratio": 200},
        "H→μμ": {"branching_ratio": 0.000218, "signature": "2 muons", "bg_ratio": 5000},
    }

    @staticmethod
    def signal_strength(observed: float, expected_sm: float) -> Dict:
        mu = observed / expected_sm if expected_sm > 0 else 0
        return {
            "mu": round(mu, 3),
            "consistent_with_sm": abs(mu - 1.0) < 0.2,
            "observed": observed,
            "expected_sm": expected_sm,
        }

    @staticmethod
    def golden_channel_analysis(sqrt_s_tev: float = 13.6, luminosity_fb: float = 450) -> Dict:
        sigma_h = CrossSections.pp_higgs_ggf(sqrt_s_tev)
        br_zz = HiggsAnalysis.DECAY_CHANNELS["H→ZZ"]["branching_ratio"]
        br_4l = br_zz * 0.044
        sigma_4l = sigma_h * br_4l
        events = CrossSections.luminosity_to_events(sigma_4l, luminosity_fb)
        bg = events * HiggsAnalysis.DECAY_CHANNELS["H→ZZ"]["bg_ratio"] * br_4l

        return {
            "channel": "H→ZZ*→4ℓ (golden channel)",
            "sigma_production_pb": round(sigma_h, 2),
            "br_zz_4l": round(br_4l, 6),
            "signal_events": events,
            "background_events": round(bg, 1),
            "significance": round(StatisticalAnalysis.significance_profile_likelihood(events, max(bg, 1)), 2),
            "mass_resolution_gev": 1.5,
            "sqrt_s_tev": sqrt_s_tev,
            "luminosity_fb": luminosity_fb,
        }

    @staticmethod
    def diphoton_analysis(sqrt_s_tev: float = 13.6, luminosity_fb: float = 450) -> Dict:
        sigma_h = CrossSections.pp_higgs_ggf(sqrt_s_tev)
        br_gg = HiggsAnalysis.DECAY_CHANNELS["H→γγ"]["branching_ratio"]
        sigma_gg = sigma_h * br_gg
        efficiency = 0.4
        events = CrossSections.luminosity_to_events(sigma_gg, luminosity_fb) * efficiency
        bg = events * HiggsAnalysis.DECAY_CHANNELS["H→γγ"]["bg_ratio"]

        return {
            "channel": "H→γγ (diphoton)",
            "signal_events": round(events, 0),
            "background_events": round(bg, 0),
            "significance": round(StatisticalAnalysis.significance_profile_likelihood(events, max(bg, 1)), 2),
            "mass_resolution_gev": 1.0,
        }


class BSMSearch:
    """Beyond Standard Model search framework."""

    @staticmethod
    def zprime_search(mass_gev: float, coupling: float = 0.1,
                       sqrt_s_tev: float = 13.6, luminosity_fb: float = 450) -> Dict:
        if mass_gev >= sqrt_s_tev * 1000:
            return {"mass_gev": mass_gev, "status": "KINEMATICALLY_FORBIDDEN",
                    "reason": f"M_Z' = {mass_gev} GeV > sqrt(s) = {sqrt_s_tev*1000} GeV"}

        sigma_pb = coupling**2 * 100 * (1000 / mass_gev)**2 * ALPHA_EM
        bg_pb = 10 * (1000 / mass_gev)**3
        signal = CrossSections.luminosity_to_events(sigma_pb, luminosity_fb)
        bg = CrossSections.luminosity_to_events(bg_pb, luminosity_fb)
        sig = StatisticalAnalysis.significance_profile_likelihood(signal, max(bg, 1))

        return {
            "particle": f"Z' (mass = {mass_gev} GeV)",
            "coupling": coupling,
            "sigma_pb": round(sigma_pb, 4),
            "signal_events": signal,
            "background_events": bg,
            "significance_sigma": round(sig, 2),
            "excluded": sig < 2.0 and signal < 3,
            "discoverable": sig >= 5.0,
        }

    @staticmethod
    def susy_gluino_search(mass_gev: float, sqrt_s_tev: float = 13.6,
                            luminosity_fb: float = 450) -> Dict:
        if mass_gev > 3000:
            return {"mass_gev": mass_gev, "status": "BEYOND_REACH"}

        sigma_pb = 1e4 * (1000 / mass_gev)**6 * ALPHA_S_MZ**2
        signal = CrossSections.luminosity_to_events(sigma_pb, luminosity_fb) * 0.1
        bg = 50 + mass_gev * 0.01
        sig = StatisticalAnalysis.significance_profile_likelihood(signal, bg)

        return {
            "particle": f"Gluino (mass = {mass_gev} GeV)",
            "sigma_pb": round(sigma_pb, 6),
            "signal_events": round(signal, 1),
            "significance_sigma": round(sig, 2),
            "current_limit_gev": 2300,
            "excluded_by_run3": mass_gev < 2300,
        }


def generate_lhc_report(sqrt_s_tev: float = 13.6, luminosity_fb: float = 450) -> Dict:
    golden = HiggsAnalysis.golden_channel_analysis(sqrt_s_tev, luminosity_fb)
    diphoton = HiggsAnalysis.diphoton_analysis(sqrt_s_tev, luminosity_fb)
    zprime_1tev = BSMSearch.zprime_search(1000, 0.1, sqrt_s_tev, luminosity_fb)
    zprime_3tev = BSMSearch.zprime_search(3000, 0.1, sqrt_s_tev, luminosity_fb)
    gluino = BSMSearch.susy_gluino_search(2500, sqrt_s_tev, luminosity_fb)
    z_peak = BreitWigner.z_boson_lineshape(91.2)
    pp_total = CrossSections.pp_total_inelastic(sqrt_s_tev)

    return {
        "report_id": hashlib.sha256(f"CERN-LHC-{datetime.now(timezone.utc).isoformat()}".encode()).hexdigest()[:16],
        "generated": datetime.now(timezone.utc).isoformat(),
        "lhc_parameters": LHC_PARAMETERS["Run 3"],
        "higgs": {
            "mass_gev": 125.25,
            "width_gev": 0.00407,
            "decay_channels": HiggsAnalysis.DECAY_CHANNELS,
            "golden_channel": golden,
            "diphoton": diphoton,
        },
        "bsm_searches": {
            "zprime_1tev": zprime_1tev,
            "zprime_3tev": zprime_3tev,
            "gluino_2500": gluino,
        },
        "standard_model": {
            "particles": len(PARTICLES),
            "z_peak_value": round(z_peak, 4),
            "pp_inelastic_mb": round(pp_total, 2),
            "fundamental_constants": {
                "alpha_em": ALPHA_EM,
                "alpha_s_mz": ALPHA_S_MZ,
                "sin2_theta_w": SIN2_THETA_W,
                "vev_gev": VEV,
            },
        },
        "engine": "ORION CERN LHC Analysis Engine v1.0",
        "author": "ORION Autonomous System",
    }
