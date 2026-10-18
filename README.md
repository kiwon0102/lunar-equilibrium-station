# Lunar Equilibrium Station — supporting computation

Reproducible calculations behind a Phase 1 concept submission to NASA's
[Lunar Grounding Challenge](https://work.crowdplat.com/challenge/lunar-grounding-challenge):
a lander-mounted station that discharges a suited astronaut from −100 kV to lander
potential without contact and without an arc.

Everything here is computation, not measurement. No experimental result in the
submission was obtained by the author; every quantitative claim is cited to a source,
computed from stated inputs by the code in this repository, or identified as an
estimate.

## What is here

| File | What it does |
| --- | --- |
| `bem_capacitance.py` | Boundary-element capacitance solver, with its validation harness |
| `compliance_case.py` | Every charge-dependent number in the concept, at the mandated 200 pF / 110 kV model |
| `data/compliance_data.csv` | Output of `compliance_case.py` — 68 values with section, unit and basis |

## Running it

```bash
pip install -r requirements.txt
python3 bem_capacitance.py     # validation, then the submission configurations
python3 compliance_case.py     # regenerates data/compliance_data.csv
```

Both scripts are self-contained and depend only on NumPy.

## The capacitance solver

The challenge specifies a 200 pF capacitance model. A first-order check treats the
crew member as an isolated conducting sphere, giving about 111 pF at 1 m radius, but
that ignores the regolith dielectric half-space underfoot and the lander conductor
alongside. How far wrong the approximation is was an open risk in the concept, so it
was computed.

Method of moments: conducting surfaces are divided into patches each carrying an
unknown charge, requiring equal potential on every patch gives a linear system, and
capacitance follows from the total charge at unit potential. The regolith and the
lander enter as image charges — a grounded conductor at coefficient −1, a dielectric
half-space at −(εr−1)/(εr+1).

That sign was the one error found during development. A positive dielectric
coefficient produces capacitances *below* the free-space value, which cannot happen:
adding a polarisable medium beneath a conductor can only raise its capacitance. The
sign was caught by that check, not by inspection.

### Validation

The solver is checked against two closed-form geometries before being used. Running
`bem_capacitance.py` prints:

```
Validation 1 - isolated conducting sphere, C = 4 pi eps0 a      worst error 0.210%
Validation 2 - sphere above a grounded plane, image series      worst error 0.254%
limit check at h/a = 50: 112.6 pF against free-space 111.3 pF
```

### Result

Crew geometry: a 0.32 m radius, 1.5 m tall cylinder for torso and life-support pack
with a 0.22 m radius spherical helmet, 1,988 patches.

| Configuration | Capacitance | vs the 200 pF model |
| --- | --- | --- |
| Free space | 64 pF | 0.32× |
| On regolith, εr = 3.6 | 88 pF | 0.44× |
| On regolith, εr = 4.22 | 92 pF | 0.46× |
| 3 m from lander | 104 pF | 0.52× |
| 2 m from lander | 113 pF | 0.56× |
| 1 m from lander | 147 pF | 0.73× |
| 0.5 m from lander | 329 pF | 1.64× |

At realistic standing distances the computed capacitance is below the mandated model,
so the challenge figure is conservative for this geometry. Even the worst case passes
the 30 s requirement with 1.8× margin.

## The compliance case

`compliance_case.py` recomputes every charge-dependent quantity against the model the
challenge mandates — 200 pF at a 110 kV differential — rather than the 111 pF sphere
estimate an earlier revision used. Q = 22.0 µC, stored energy 1.21 J, 733 nA average
over 30 s.

Two results from that pass changed the design:

- **The Stage C bleed resistor.** With the platform capacitance of 4,427 nF carried
  through the time constant, 10 MΩ gives 5τ = 221 s and breaks the 30 s requirement.
  The baseline moved to a standard 470 kΩ part: 5τ = 10.4 s at 10.6 µA.
- **Regolith as a return path is rejected.** It would need below 1.0 × 10¹⁰ Ω·m.
  Apollo-era measurements put regolith at 10¹³–10¹⁶ Ω·m, reaching 9.1 × 10¹⁶ Ω·m at
  the −218 °C floor via σ(T) = 3 × 10⁻¹⁸ exp(0.0237 T) from sample 15301,38.

## Sources

Apollo-era regolith conductivity: *Science* (1969), *Earth and Planetary Science
Letters* (1972), *Journal of Geophysical Research* (1974, sample 15301,38). Simulant
electrostatics: NASA KSC, NTRS 20130012059. Dust shield flight result: NTRS
20250009407. JSC-1A permittivity 3.6–4.22. Full reference list is in the submission
document.

## AI disclosure

Claude (Anthropic) was used to draft the submission documents, search the literature
and write the code in this repository, under the author's direction and review. The
author set the design direction, made the engineering decisions and verified the
results. The code is published here so the computations can be checked independently
rather than taken on trust.

## Licence

MIT. See `LICENSE`.
