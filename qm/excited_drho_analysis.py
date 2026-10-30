#!/usr/bin/env python3
"""Excited-state analogue of cube_drho.py: charge-displacement analysis and
excited-state ESP from TD-DFT difference densities (orca_plot option 6,
Δρ = ρ_exc − ρ_0, unrelaxed) of a two-fragment complex.

    python3 excited_drho_analysis.py GSDIR NA label1:cube1 [label2:cube2 ...]

GSDIR holds the ground-state cubes from make_charge_density_cubes.sh
(complex_dens.cub, complex_esp.cub, <fragA>_dens.cub, <fragB>_dens.cub given as
rubpy_dens.cub / mos_dens.cub). Fragment A = first NA atoms (acceptor side in
the CD axis), B = the rest.

Per state it prints/writes:
  * q_CT(boundary): Δq(t) at the ground-state isodensity boundary (charge-
    displacement curve, Belpassi JACS 2008); negative = electrons B -> A
  * Δq_A: Δρ integrated over the region where ρ_A > ρ_B (ground state)
  * q_gain = ∫Δρ+ and D = distance between the centroids of Δρ+ and Δρ−
    (Le Bahers JCTC 2011)
  * ΔV = potential of Δρ (FFT convolution with 1/r on a 2x-coarsened grid),
    and <label>_esp.cub = ground-state ESP + ΔV on the ground-state grid
Writes <label>_cd.dat, excited_cd.{png,pdf} and excited_summary.txt.
"""
import sys
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import fftconvolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cube_drho import read_cube, write_cube, BOHR  # noqa: E402

COLORS = ["#21b0fe", "#fe218b", "#fed700", "#2ca02c", "#7f3fbf", "#ff7f0e"]


def axes_coords(c):
    return [c["org"][i] + np.arange(c["n"][i]) * c["ax"][i][i] for i in range(3)]


def interp_to(src, dst_pts):
    f = RegularGridInterpolator(axes_coords(src), src["data"], bounds_error=False, fill_value=0.0)
    return f(dst_pts.reshape(-1, 3)).reshape(dst_pts.shape[:-1])


def grid_points(c):
    x, y, z = axes_coords(c)
    return np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1)


def delta_potential(c, d):
    """Potential (a.u.) of the electron density difference d: V = -∫ d(r')/|r-r'|."""
    n = [m - m % 2 for m in c["n"]]
    dd = d[:n[0], :n[1], :n[2]].reshape(n[0] // 2, 2, n[1] // 2, 2, n[2] // 2, 2).sum(axis=(1, 3, 5))
    h = np.array([2 * c["ax"][i][i] for i in range(3)])
    dV0 = np.prod(h / 2)
    q = dd * dV0  # charge per coarse cell
    m = dd.shape
    r = [np.arange(-(k - 1), k) * h[i] for i, k in enumerate(m)]
    X, Y, Z = np.meshgrid(*r, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    heff = np.prod(h) ** (1 / 3)
    with np.errstate(divide="ignore"):
        G = np.where(R > 0, 1.0 / R, 2.380 / heff)
    V = -fftconvolve(q, G, mode="same")
    xc = [c["org"][i] + (np.arange(m[i]) * 2 + 0.5) * c["ax"][i][i] for i in range(3)]
    return RegularGridInterpolator(xc, V, bounds_error=False, fill_value=None)


def main():
    gsdir, na = sys.argv[1], int(sys.argv[2])
    states = [a.split(":", 1) for a in sys.argv[3:]]
    gs = read_cube(f"{gsdir}/complex_dens.cub")
    gA = read_cube(f"{gsdir}/rubpy_dens.cub")
    gB = read_cube(f"{gsdir}/mos_dens.cub")
    gesp = read_cube(f"{gsdir}/complex_esp.cub")
    xyz = gs["xyz"]
    cA, cB = xyz[:na].mean(0), xyz[na:].mean(0)
    u = (cA - cB) / np.linalg.norm(cA - cB)
    gpts = grid_points(gs)
    surf = (gs["data"] > 0.0009) & (gs["data"] < 0.0011)
    surfA = surf & (gA["data"] > gB["data"])

    # ground-state isodensity boundary along the axis (same definition as cube_drho.py)
    tg = ((gpts - cB) @ u) * BOHR
    edges = np.arange(tg.min(), tg.max() + 0.05, 0.05)
    rA, _ = np.histogram(tg, edges, weights=gA["data"])
    rB, _ = np.histogram(tg, edges, weights=gB["data"])
    tc = edges[1:]
    L = np.linalg.norm(cA - cB) * BOHR
    win = (tc > 0) & (tc < L)
    tb = tc[win][np.argmin(np.abs(rA[win] - rB[win]))]

    rows, curves = [], []
    for lab, fn in states:
        c = read_cube(fn)
        d = c["data"]
        dV = abs(np.linalg.det(c["ax"]))
        pts = grid_points(c)
        t = ((pts - cB) @ u) * BOHR
        e2 = np.arange(t.min(), t.max() + 0.05, 0.05)
        dq = np.cumsum(np.histogram(t, e2, weights=d * dV)[0])
        t2 = e2[1:]
        qb = np.interp(tb, t2, dq)
        inA = interp_to(gA, pts) > interp_to(gB, pts)
        dqA = d[inA].sum() * dV
        pos, neg = np.clip(d, 0, None), np.clip(-d, 0, None)
        qg, ql = pos.sum() * dV, neg.sum() * dV
        cpos = (pts * pos[..., None]).reshape(-1, 3).sum(0) / pos.sum()
        cneg = (pts * neg[..., None]).reshape(-1, 3).sum(0) / neg.sum()
        D = np.linalg.norm(cpos - cneg) * BOHR
        dvf = delta_potential(c, d)
        dVg = dvf(gpts.reshape(-1, 3)).reshape(gs["data"].shape)
        esp_exc = gesp["data"] + dVg
        write_cube(f"{lab}_esp.cub", gs, esp_exc, f"{lab} ESP = ground-state ESP + potential of TD-DFT difference density")
        write_cube(f"{lab}_dV.cub", gs, dVg, f"{lab} potential of TD-DFT difference density")
        k = 627.5095
        sA, sB = esp_exc[surfA] * k, esp_exc[surf & ~surfA] * k
        dA, dB = dVg[surfA] * k, dVg[surf & ~surfA] * k
        np.savetxt(f"{lab}_cd.dat", np.c_[t2, dq], header="t/A (0 = centroid B)  Delta_q(t)/e", fmt="%.4f")
        curves.append((lab, t2, dq))
        rows.append((lab, int(np.sum(d)), dq[-1], qb, dqA, qg, ql, D,
                     sA.mean(), sB.mean(), dA.mean(), dB.mean(), np.percentile(np.abs(esp_exc[surf]), 98)))
        print(f"{lab}: int Δρ {d.sum()*dV:+.4f}  q_CT(boundary) {qb:+.3f}  Δq_A {dqA:+.3f}  "
              f"q_gain {qg:.3f}  D {D:.2f} A  ESP_A {sA.mean():+.1f} ESP_B {sB.mean():+.1f}  "
              f"ΔV_A {dA.mean():+.1f} ΔV_B {dB.mean():+.1f} kcal/mol·e", flush=True)

    with open("excited_summary.txt", "w") as f:
        f.write(f"# boundary t_b = {tb:.2f} A from centroid B; centroid distance {L:.2f} A\n")
        f.write("# state  int_drho  q_CT(boundary,e; <0 = B->A)  dq_A(e)  q_gain  q_loss  D(A)  "
                "<ESP>_A  <ESP>_B  <dV>_A  <dV>_B (kcal/mol/e)  |ESP|98\n")
        for r in rows:
            f.write(f"{r[0]:10s} {r[2]:+.4f} {r[3]:+.3f} {r[4]:+.3f} {r[5]:.3f} {r[6]:.3f} {r[7]:.2f} "
                    f"{r[8]:+.1f} {r[9]:+.1f} {r[10]:+.1f} {r[11]:+.1f} {r[12]:.4f}\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    ax.axhline(0, color="0.6", lw=0.7)
    ax.axvline(tb, color="0.4", ls="--", lw=0.9)
    for (lab, t2, dq), col in zip(curves, COLORS):
        ax.plot(t2, dq, lw=1.4, color=col, label=lab.replace("_", " "))
    tA = (xyz[:na] - cB) @ u * BOHR
    tB = (xyz[na:] - cB) @ u * BOHR
    y0 = ax.get_ylim()[0]
    ax.plot(tB, np.full_like(tB, y0), "|", color="0.55", ms=5)
    ax.plot(tA, np.full_like(tA, y0), "|", color="0.25", ms=5)
    ax.set_xlabel("position along Mo$_3$S$_{13}$$\\rightarrow$Ru(bpy)$_3$ axis (Å)")
    ax.set_ylabel("$\\Delta q(t)$ (e)")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"excited_cd.{ext}", dpi=300)


if __name__ == "__main__":
    main()
