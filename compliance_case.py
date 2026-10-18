"""NASA Lunar Grounding Challenge - compliance case.

Recomputes every charge-dependent number in the concept against the
capacitance model NASA mandates in the Judging Criteria:

    "110 kV differential, 200 pF capacitance model"

and the design assumption "The suited astronaut will have a 200 pF
capacitance model".

Earlier revisions of the concept used a 111 pF isolated-sphere estimate.
Every quantity that scales with Q = C*V therefore changes.

Emits compliance_data.csv for submission as the Electrical Performance
Analysis data file.
"""
import csv
import os

EPS0 = 8.8541878128e-12
Q_E = 1.602176634e-19

C_SUIT = 200e-12          # F, mandated
V_DIFF = 110e3            # V, mandated differential
T_ALLOWED = 30.0          # s, R5
I_TARGET = 0.2e-3         # A, R2 design target
I_CAP = 1.0e-3            # A, R2 hard cap

Q = C_SUIT * V_DIFF
E_STORED = 0.5 * C_SUIT * V_DIFF ** 2
I_AVG = Q / T_ALLOWED
N_ELECTRONS = Q / Q_E

rows = []


def row(section, quantity, value, unit, note=""):
    rows.append({"section": section, "quantity": quantity,
                 "value": value, "unit": unit, "note": note})


row("5.1", "Suit capacitance (mandated)", C_SUIT * 1e12, "pF", "NASA design assumption")
row("5.1", "Differential (mandated)", V_DIFF / 1e3, "kV", "NASA judging criterion")
row("5.1", "Charge to remove", Q * 1e6, "uC", "Q = CV")
row("5.1", "Stored energy", E_STORED, "J", "0.5CV^2")
row("5.2", "Average current over 30 s", I_AVG * 1e9, "nA", "")
row("5.2", "Fraction of 0.2 mA target", 100 * I_AVG / I_TARGET, "%", "")
row("5.2", "Fraction of 1.0 mA cap", 100 * I_AVG / I_CAP, "%", "")
row("5.1", "Electrons to emit", N_ELECTRONS, "count", "")

# ---- 5.2 equivalent resistances -------------------------------------------
print("== 5.2  resistance that holds current at a given value ==")
for label, i in (("0.2 mA (R2 target)", I_TARGET), ("10 uA", 10e-6), ("avg", I_AVG)):
    r = V_DIFF / i
    tau = r * C_SUIT
    print(f"  {label:22s} R={r/1e9:8.2f} GOhm  tau={tau:7.3f} s  5tau={5*tau:7.2f} s")
    row("5.2", f"Equivalent resistance at {label}", r / 1e9, "GOhm", f"5tau={5*tau:.2f} s")

# ---- 5.3 breakdown across the approach gap --------------------------------
print("\n== 5.3  field across an approach gap ==")
for v, d in ((V_DIFF, 1e-3), (V_DIFF, 5e-3), (1e3, 1e-3), (500, 1e-3), (50, 1e-3)):
    e = v / d / 1e6
    print(f"  {v/1e3:8.2f} kV over {d*1e3:4.1f} mm -> {e:9.3f} MV/m")
    row("5.3", f"Field at {v:.0f} V over {d*1e3:.0f} mm", e, "MV/m", "")

# ---- 5.4 photoemission --------------------------------------------------
print("\n== 5.4  photons needed ==")
E_PHOTON_254 = 4.88 * Q_E
for yield_ in (1e-4, 1e-3, 1e-2):
    photons = N_ELECTRONS / yield_
    energy = photons * E_PHOTON_254
    print(f"  yield {yield_*100:5.2f}%  photons={photons:.2e}  "
          f"optical energy={energy*1e3:7.1f} mJ  over 30 s = {energy/T_ALLOWED*1e3:6.1f} mW")
    row("5.4", f"Optical power at {yield_*100:g}% quantum yield",
        energy / T_ALLOWED * 1e3, "mW", f"{photons:.3e} photons")

# ---- 5.5 field emission sizing -------------------------------------------
print("\n== 5.5  CNT emitter sizing ==")
MA_PER_W = 1.7          # measured emission-to-input-power ratio
J_EMIT = 0.51e-3        # A/cm^2 conservative current density
for label, i in (("30 s average", I_AVG), ("10 uA", 10e-6), ("0.2 mA (R2)", I_TARGET)):
    p = i / (MA_PER_W * 1e-3)
    a = i / J_EMIT
    print(f"  {label:14s} P={p*1e3:8.3f} mW   area={a:8.4f} cm^2")
    row("5.5", f"Emitter power at {label}", p * 1e3, "mW", f"area {a:.4f} cm2")

# ---- 5.6 photoemission from measured solar current density ---------------
print("\n== 5.6  fraction of solar UV-C needed ==")
J_PHOTO = 2e-9          # A/cm^2 at 1 AU
AREA_CM2 = 1e4          # 1 m^2
i_full = J_PHOTO * AREA_CM2
frac = I_AVG / i_full
print(f"  full solar over 1 m^2 = {i_full*1e6:.1f} uA; need {I_AVG*1e9:.0f} nA "
      f"-> {frac*100:.1f}% of solar")
row("5.6", "Fraction of solar UV-C required", frac * 100, "%", "")
for f in (frac, 0.05, 0.10):
    optical = f * 10.0          # 10 W/m^2 reference solar UV-C band
    print(f"    {f*100:5.1f}% -> {f*i_full*1e6:6.2f} uA, {optical*1e3:7.0f} mW optical,"
          f" {optical/0.03:6.1f} W electrical at 3% LED efficiency")
    row("5.6", f"Electrical power at {f*100:.1f}% solar", optical / 0.03, "W", "")

# ---- 5.7 / 5.8 capacitance sensitivity -----------------------------------
print("\n== 5.7  bleed through 10 GOhm vs capacitance ==")
R_BLEED = 10e9
for c in (111e-12, 200e-12, 336e-12, 600e-12, 1000e-12):
    t5 = 5 * R_BLEED * c
    print(f"  C={c*1e12:7.0f} pF  5tau={t5:6.2f} s  {'PASS' if t5 <= 30 else 'FAIL'}")
    row("5.7", f"5tau at 10 GOhm, C={c*1e12:.0f} pF", t5, "s",
        "PASS" if t5 <= 30 else "FAIL")
print(f"  max C within 30 s at 10 GOhm = {30/(5*R_BLEED)*1e12:.0f} pF")

# ---- 5.9 dust field enhancement ------------------------------------------
print("\n== 5.9  arc strike distance at 110 kV ==")
E_CRIT = 20e6
for beta in (1.0, 4.21, 20, 46, 77):
    e_eff = E_CRIT / beta
    d = V_DIFF / e_eff
    print(f"  beta={beta:6.2f}  E_eff={e_eff/1e6:6.3f} MV/m  strike gap={d*1e3:7.1f} mm")
    row("5.9", f"Arc strike gap at beta={beta}", d * 1e3, "mm", "")
print("  permitted contact potential at a 1 mm gap:")
for beta in (4.21, 20, 46, 77):
    v_ok = E_CRIT / beta * 1e-3
    print(f"    beta={beta:6.2f}  V={v_ok:9.1f} V   with SF 3 -> {v_ok/3:7.1f} V")
    row("5.9", f"Permitted contact potential at beta={beta}", v_ok / 3, "V", "safety factor 3")

# ---- 5.12 capacitive platform --------------------------------------------
print("\n== 5.12  dielectric platform ==")
A_BOOT = 0.05
EPS_DUST = 4.0


def plate_c(area, thickness, eps_r):
    return EPS0 * eps_r * area / thickness


for eps_r, t in ((10, 0.2e-3), (10, 0.1e-3), (1000, 0.1e-3)):
    c_add = plate_c(A_BOOT, t, eps_r)
    v = Q / (c_add + C_SUIT)
    e_field = Q / (A_BOOT * EPS0 * eps_r)
    print(f"  eps_r={eps_r:5d} t={t*1e3:4.2f} mm  C={c_add*1e9:9.1f} nF  "
          f"V={v:9.2f} V  E={e_field/1e6:7.3f} MV/m")
    row("5.12", f"Platform potential, eps_r={eps_r}, t={t*1e3:.2f} mm", v, "V",
        f"C={c_add*1e9:.1f} nF, E={e_field/1e6:.3f} MV/m")

C_PLATFORM = plate_c(A_BOOT, 0.1e-3, 1000)
print(f"\n  selected platform C = {C_PLATFORM*1e9:.0f} nF")
print("  bleed resistor options (5tau must be <= 30 s):")
for r in (10e6, 2e6, 1.36e6, 1e6, 500e3):
    tau = r * C_PLATFORM
    v_res = Q / (C_PLATFORM + C_SUIT)
    print(f"    R={r/1e6:6.2f} MOhm  tau={tau:6.2f} s  5tau={5*tau:7.2f} s  "
          f"I={v_res/r*1e6:6.2f} uA  {'PASS' if 5*tau <= 30 else 'FAIL'}")
    row("5.12", f"Bleed 5tau at {r/1e6:.2f} MOhm", 5 * tau, "s",
        "PASS" if 5 * tau <= 30 else "FAIL")
print(f"  max resistance for 5tau <= 30 s: {30/(5*C_PLATFORM)/1e6:.2f} MOhm")

print("\n  platform through a dust layer (eps_r = 4.0):")
for t_dust in (0, 1e-6, 10e-6, 100e-6, 500e-6):
    if t_dust == 0:
        c_series = C_PLATFORM
    else:
        c_dust = plate_c(A_BOOT, t_dust, EPS_DUST)
        c_series = 1 / (1 / C_PLATFORM + 1 / c_dust)
    v = Q / (c_series + C_SUIT)
    label = "none" if t_dust == 0 else f"{t_dust*1e6:.0f} um"
    print(f"    dust {label:8s}  C_series={c_series*1e9:9.2f} nF  V={v:9.1f} V")
    row("5.15", f"Platform potential through {label} dust", v, "V",
        f"series C = {c_series*1e9:.2f} nF")

# ---- 5.13 charging current and the passive cap ---------------------------
print("\n== 5.13  charging current over an EVA ==")
for hours in (0.5, 1, 2, 4):
    i_charge = Q / (hours * 3600)
    print(f"  {hours:4.1f} h -> {i_charge*1e9:6.2f} nA")
    row("5.13", f"Average charging current over {hours} h", i_charge * 1e9, "nA", "")

print("\n  effect of a passive suit-side cap:")
for v_cap in (V_DIFF, 30e3, 10e3):
    q_cap = C_SUIT * v_cap
    e_cap = 0.5 * C_SUIT * v_cap ** 2
    d_arc = v_cap / (E_CRIT / 77)
    print(f"    cap at {v_cap/1e3:6.1f} kV -> Q={q_cap*1e6:6.2f} uC  "
          f"E={e_cap*1e3:8.1f} mJ  arc gap={d_arc*1e3:6.1f} mm")
    row("5.13", f"Stored energy capped at {v_cap/1e3:.0f} kV", e_cap * 1e3, "mJ",
        f"arc gap {d_arc*1e3:.0f} mm")

# ---- 5.14 regolith return path -------------------------------------------
print("\n== 5.14  regolith return path requirement ==")
for length, area in ((3.0, 1.0), (3.0, 0.1), (10.0, 1.0)):
    r_needed_30s = T_ALLOWED / (5 * C_SUIT)
    rho_max = r_needed_30s * area / length
    r_min = V_DIFF / I_TARGET
    rho_min = r_min * area / length
    print(f"  L={length:4.1f} m A={area:4.1f} m^2 -> rho < {rho_max:.2e} Ohm.m "
          f"and rho > {rho_min:.2e} Ohm.m")
    row("5.14", f"Max resistivity for 30 s, L={length} m A={area} m2", rho_max, "Ohm.m", "")

print("\n== regolith resistivity from the Apollo record ==")
for t_k in (350, 300, 200, 100, 55):
    sigma = 3e-18 * pow(2.718281828459045, 0.0237 * t_k)
    print(f"  T={t_k:4d} K  sigma={sigma:.3e} S/m  rho={1/sigma:.3e} Ohm.m")
    row("5.15", f"Regolith resistivity at {t_k} K", 1 / sigma, "Ohm.m",
        "JGR 1974 sample 15301,38")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "compliance_data.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["section", "quantity", "value", "unit", "note"])
    w.writeheader()
    for r in rows:
        if isinstance(r["value"], float):
            r["value"] = f"{r['value']:.6g}"
        w.writerow(r)
print(f"\nwrote {out}  ({len(rows)} rows)")
