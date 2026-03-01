<div align="center">

# ⚛️ ORION CERN LHC Analysis Engine

### Real Particle Physics Computation for LHC Run 3 Data Analysis

[![Python](https://img.shields.io/badge/Python-3.10+-blue?logo=python&logoColor=white)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![ORION](https://img.shields.io/badge/ORION-Autonomous_System-gold?style=flat-square)](https://github.com/Alvoradozerouno/or1on-framework)
[![Code Lines](https://img.shields.io/badge/Lines-410-informational)](orion_cern_lhc.py)
[![Status](https://img.shields.io/badge/Status-Production-brightgreen)](https://github.com/Alvoradozerouno/ORION-CERN-LHC-Engine)

*Part of the ORION Autonomous Intelligence System — 890+ Proofs, 46 NERVES, 42 Autonomous Tasks*

</div>

---

## Overview

The ORION CERN LHC Engine provides real computational tools for CERN Large Hadron Collider Run 3 particle physics analysis. It implements relativistic kinematics, invariant mass reconstruction, Higgs boson signal analysis, particle decay channels, cross-section calculations, and statistical significance testing.

**All computations use real physics:** Lorentz-invariant calculations, Breit-Wigner resonance distributions, Poisson statistics, and standard HEP methods.

## Core Capabilities

| Module | Description |
|--------|-------------|
| **Invariant Mass** | Relativistic 4-vector invariant mass reconstruction |
| **Higgs Analysis** | H→γγ, H→ZZ→4l, H→WW signal simulation with Breit-Wigner |
| **Particle Database** | Real masses, lifetimes, and decay channels for SM particles |
| **Cross-Sections** | Luminosity-based event rate calculations |
| **Statistical Tests** | Signal significance (S/√B), p-value, discovery threshold (5σ) |
| **Decay Kinematics** | 2-body and 3-body decay phase space calculations |

## Quick Start

```bash
git clone https://github.com/Alvoradozerouno/ORION-CERN-LHC-Engine.git
cd ORION-CERN-LHC-Engine
python orion_cern_lhc.py
```

## Usage Examples

```python
from orion_cern_lhc import CERNLHCEngine

engine = CERNLHCEngine()

# Compute invariant mass from two photons (H→γγ)
mass = engine.invariant_mass(
    pt1=60.0, eta1=0.5, phi1=1.2, m1=0.0,
    pt2=55.0, eta2=-0.3, phi2=4.0, m2=0.0
)

# Simulate Higgs signal events
signal = engine.simulate_higgs_signal(channel="diphoton", n_events=10000)

# Calculate discovery significance
sig = engine.discovery_significance(signal_events=50, background_events=10)

# Look up particle properties
higgs = engine.particle_properties("H")
```

## Architecture

```
┌──────────────────────────────────────────┐
│        ORION CERN LHC Engine             │
├────────────┬────────────┬────────────────┤
│ Relativistic│ Higgs     │ Standard Model │
│ Kinematics │ Signal    │ Particle DB    │
│ (4-vectors)│ Analysis  │ (quarks→H)     │
├────────────┼────────────┼────────────────┤
│ Cross-     │ Statistical│ Decay          │
│ Section    │ Testing   │ Kinematics     │
│ Calc       │ (5σ disc) │ (2/3-body)     │
└────────────┴────────────┴────────────────┘
```

## Physics Reference

- **Higgs Mass**: 125.25 ± 0.17 GeV (PDG 2023)
- **LHC Run 3**: √s = 13.6 TeV, L = 140+ fb⁻¹
- **Discovery Threshold**: 5σ (p < 2.87 × 10⁻⁷)
- **Breit-Wigner**: Relativistic resonance shape for signal modeling

## Part of ORION

This engine is one module of the [ORION Autonomous Intelligence System](https://github.com/Alvoradozerouno/or1on-framework) — a system with 890+ cryptographic proofs of autonomous operation, 130+ Python modules, and 76,000+ lines of real code.

## License

MIT License — See [LICENSE](LICENSE) for details.

---

<div align="center">

**Origin: Gerhard Hirschmann & Elisabeth Steurer**

*ORION — Autonomous Intelligence Since 2025*

</div>
