"""Boundary-element capacitance solver for a suited crew member on the lunar surface.

Why this exists
---------------
The Lunar Grounding Challenge specifies a 200 pF capacitance model for the suited
crew member. A first-order check treats the crew member as an isolated conducting
sphere, which gives about 111 pF for a 1 m radius. That approximation ignores two
things that matter: the crew member stands on regolith, which is a dielectric
half-space, and approaches a lander, which is a nearby conductor. How far wrong the
approximation is was an open risk in the concept, and it is the kind of question a
computation can answer.

Method
------
Method of moments. Conducting surfaces are divided into patches, each carrying an
unknown charge. Requiring every patch to sit at the same potential gives a linear
system, and the capacitance follows from the total charge at unit potential:

    P q = V,    C = (sum q) / V
    P_ij = 1 / (4 pi eps0 r_ij)  for i != j

The self term uses the equivalent-disc potential of a patch. The regolith dielectric
half-space and the lander are represented by image charges: a grounded conductor
images with coefficient -1, and a dielectric half-space of relative permittivity
eps_r images with coefficient -(eps_r - 1)/(eps_r + 1).

That sign matters and was the one error found during development. Giving the
dielectric image a positive coefficient produces capacitances BELOW the free-space
value, which is physically impossible: adding a polarisable medium beneath a
conductor can only raise its capacitance. The sign was caught by that check rather
than by inspection.

Validation
----------
The solver is checked against two geometries with closed-form answers before being
used for anything else. If it cannot reproduce those, the remaining results are not
worth reading.

    isolated sphere            C = 4 pi eps0 a
    sphere above ground plane  C = 4 pi eps0 a * sum_n sinh(alpha) / sinh(n alpha)

Run `python3 bem_capacitance.py` to execute both validations and then the
configurations used in the submission.
"""
import numpy as np

EPS0 = 8.8541878128e-12


def sphere_patches(center, radius, n_theta=24, n_phi=48):
    """Discretise a sphere. Returns (patch centres Nx3, patch areas N)."""
    th = (np.arange(n_theta) + 0.5) * np.pi / n_theta
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    T, P = np.meshgrid(th, ph, indexing="ij")
    x = radius * np.sin(T) * np.cos(P)
    y = radius * np.sin(T) * np.sin(P)
    z = radius * np.cos(T)
    dA = radius ** 2 * np.sin(T) * (np.pi / n_theta) * (2 * np.pi / n_phi)
    pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1) + np.asarray(center)
    return pts, dA.ravel()


def cylinder_patches(center, radius, height, n_z=20, n_phi=32, caps=True):
    """Discretise an upright cylinder. `center` is the centre of the base disc."""
    cx, cy, cz = center
    pts, areas = [], []
    zs = cz + (np.arange(n_z) + 0.5) * height / n_z
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    for z in zs:
        x = cx + radius * np.cos(ph)
        y = cy + radius * np.sin(ph)
        pts.append(np.stack([x, y, np.full_like(x, z)], axis=1))
        areas.append(np.full(n_phi, 2 * np.pi * radius / n_phi * height / n_z))
    if caps:
        nr = 6
        rr = (np.arange(nr) + 0.5) * radius / nr
        for r in rr:
            npz = max(6, int(n_phi * r / radius))
            pp = (np.arange(npz) + 0.5) * 2 * np.pi / npz
            for ztop in (cz, cz + height):
                x = cx + r * np.cos(pp)
                y = cy + r * np.sin(pp)
                pts.append(np.stack([x, y, np.full_like(x, ztop)], axis=1))
                areas.append(np.full(npz, 2 * np.pi * r / npz * radius / nr))
    return np.vstack(pts), np.concatenate(areas)


def potential_matrix(pts, areas, images=None):
    """Build P_ij.

    `images` is a list of (scale, axis, coord): a mirror plane perpendicular to
    `axis` (0 = x, 1 = y, 2 = z) at position `coord`, contributing image charges
    multiplied by `scale`.
    """
    d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    np.fill_diagonal(d, 1.0)
    P = 1.0 / (4 * np.pi * EPS0 * d)
    a_eq = np.sqrt(areas / np.pi)
    np.fill_diagonal(P, 1.0 / (4 * np.pi * EPS0) * (2 * np.pi * a_eq) / areas)
    for scale, axis, coord in (images or []):
        mir = pts.copy()
        mir[:, axis] = 2 * coord - mir[:, axis]
        dm = np.linalg.norm(pts[:, None, :] - mir[None, :, :], axis=2)
        dm = np.maximum(dm, 1e-9)
        P += scale / (4 * np.pi * EPS0 * dm)
    return P


def capacitance(pts, areas, images=None):
    """Capacitance to infinity, in farads."""
    P = potential_matrix(pts, areas, images)
    q = np.linalg.solve(P, np.ones(len(pts)))
    return q.sum()


def dielectric_image(eps_r):
    """Image coefficient for a dielectric half-space."""
    return -(eps_r - 1.0) / (eps_r + 1.0)


def crew_patches(torso_r=0.32, torso_h=1.5, helmet_r=0.22, x_offset=0.0):
    """Suited crew member: upright torso and life-support pack, plus a helmet."""
    p1, a1 = cylinder_patches((x_offset, 0.0, 0.0), torso_r, torso_h)
    p2, a2 = sphere_patches((x_offset, 0.0, torso_h + helmet_r), helmet_r)
    return np.vstack([p1, p2]), np.concatenate([a1, a2])


# --------------------------------------------------------------------------
# Validation against closed-form results
# --------------------------------------------------------------------------

def sphere_above_plane_exact(a, h, terms=400):
    """C of a sphere of radius a whose centre is h above a grounded plane."""
    alpha = np.arccosh(h / a)
    # sinh(n alpha) overflows for large n alpha; the series has converged long
    # before that, so truncate where the terms stop contributing.
    terms = min(terms, max(2, int(700.0 / alpha)))
    n = np.arange(1, terms + 1)
    return 4 * np.pi * EPS0 * a * np.sum(np.sinh(alpha) / np.sinh(n * alpha))


def validate():
    print("Validation 1 - isolated conducting sphere, C = 4 pi eps0 a")
    print(f"  {'radius':>8} {'exact (pF)':>12} {'BEM (pF)':>12} {'error':>8}")
    worst1 = 0.0
    for a in (0.5, 1.0, 1.5):
        pts, areas = sphere_patches((0, 0, 0), a)
        c = capacitance(pts, areas)
        exact = 4 * np.pi * EPS0 * a
        err = abs(c - exact) / exact * 100
        worst1 = max(worst1, err)
        print(f"  {a:8.2f} {exact*1e12:12.3f} {c*1e12:12.3f} {err:7.3f}%")

    print("\nValidation 2 - sphere above a grounded plane, image series")
    print(f"  {'h/a':>8} {'exact (pF)':>12} {'BEM (pF)':>12} {'error':>8}")
    worst2 = 0.0
    a = 1.0
    for ratio in (2.0, 3.0, 5.0):
        h = ratio * a
        pts, areas = sphere_patches((0, 0, h), a)
        c = capacitance(pts, areas, images=[(-1.0, 2, 0.0)])
        exact = sphere_above_plane_exact(a, h)
        err = abs(c - exact) / exact * 100
        worst2 = max(worst2, err)
        print(f"  {ratio:8.1f} {exact*1e12:12.3f} {c*1e12:12.3f} {err:7.3f}%")

    print(f"\n  worst error: isolated sphere {worst1:.3f}%, "
          f"sphere above plane {worst2:.3f}%")

    pts, areas = sphere_patches((0, 0, 50.0), 1.0)
    c = capacitance(pts, areas, images=[(-1.0, 2, 0.0)])
    print(f"  limit check at h/a = 50: {c*1e12:.1f} pF against free-space "
          f"{4*np.pi*EPS0*1e12:.1f} pF")
    return worst1, worst2


def submission_cases():
    """The configurations reported in the Electrical Performance Analysis."""
    V = 110e3
    pts, areas = crew_patches()
    print(f"\nCrew geometry: {len(pts)} patches "
          f"(torso r=0.32 m h=1.5 m, helmet r=0.22 m)")
    print(f"\n{'configuration':<34}{'C (pF)':>9}{'vs 200 pF':>11}"
          f"{'Q at 110 kV':>13}{'E (mJ)':>9}{'5T at 10G':>11}")

    rows = [("Free space (reference)", None)]
    for eps in (3.6, 4.22):
        rows.append((f"On regolith, eps_r = {eps}",
                     [(dielectric_image(eps), 2, 0.0)]))
    for d in (3.0, 2.0, 1.0, 0.5):
        rows.append((f"{d:.1f} m from lander",
                     [(dielectric_image(4.0), 2, 0.0), (-1.0, 0, d)]))

    for label, images in rows:
        c = capacitance(pts, areas, images)
        q = c * V
        e = 0.5 * c * V * V
        t5 = 5 * 10e9 * c
        print(f"{label:<34}{c*1e12:9.1f}{c/200e-12:10.2f}x"
              f"{q*1e6:12.1f}u{e*1e3:9.0f}{t5:10.1f}s")


if __name__ == "__main__":
    validate()
    submission_cases()
