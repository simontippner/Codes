#!/usr/bin/env python3
"""Charge-density difference drho = rho_AB - rho_A - rho_B from three cubes on
the same grid (make_charge_density_cubes.sh), plus the charge-displacement
curve (Belpassi et al., JACS 2008, 130, 1048) along the axis joining the two
fragment centroids.

    python3 cube_drho.py complex.cub fragA.cub fragB.cub NA [prefix] [nameA] [nameB]

NA = number of atoms of fragment A (A = first NA atoms of the complex, B = rest).
Writes <prefix>.cub (drho), <prefix>_cd.dat and <prefix>_cd.{png,pdf}.

Delta q(t) = integral of drho over all points with axis coordinate < t, with t
running from fragment B (left) to fragment A (right). Delta q > 0 at the
boundary means electrons moved from A to B.
"""
import sys
import numpy as np

BOHR = 0.529177210903


def read_cube(fn):
    with open(fn) as f:
        lines = f.readlines()
    nat, *org = lines[2].split()
    nat = abs(int(nat))
    org = np.array(org[:3], float)
    n, ax = [], []
    for l in lines[3:6]:
        p = l.split()
        n.append(int(p[0]))
        ax.append(np.array(p[1:4], float))
    atoms = [l.split() for l in lines[6:6 + nat]]
    Z = np.array([int(a[0]) for a in atoms])
    xyz = np.array([a[2:5] for a in atoms], float)
    data = np.array(" ".join(lines[6 + nat:]).split(), float).reshape(n)
    return dict(head=lines[:6 + nat], org=org, n=n, ax=np.array(ax), Z=Z, xyz=xyz, data=data)


def write_cube(fn, c, data, title):
    with open(fn, "w") as f:
        f.write(title + "\n drho = rho_AB - rho_A - rho_B\n")
        f.writelines(c["head"][2:])
        flat = data.reshape(c["n"][0] * c["n"][1], c["n"][2])
        for row in flat:
            for i in range(0, len(row), 6):
                f.write("".join(f"{v:13.5E}" for v in row[i:i + 6]) + "\n")


def main():
    fab, fa, fb, na = sys.argv[1:5]
    na = int(na)
    pre = sys.argv[5] if len(sys.argv) > 5 else "drho"
    nameA = sys.argv[6] if len(sys.argv) > 6 else "A"
    nameB = sys.argv[7] if len(sys.argv) > 7 else "B"
    ab, a, b = read_cube(fab), read_cube(fa), read_cube(fb)
    for c in (a, b):
        assert c["n"] == ab["n"] and np.allclose(c["org"], ab["org"]) and np.allclose(c["ax"], ab["ax"]), "grids differ"
    d = ab["data"] - a["data"] - b["data"]
    dV = abs(np.linalg.det(ab["ax"]))
    write_cube(pre + ".cub", ab, d, "charge density difference")
    print(f"int rho_AB {ab['data'].sum()*dV:.3f}  rho_A {a['data'].sum()*dV:.3f}  rho_B {b['data'].sum()*dV:.3f}")
    print(f"int drho {d.sum()*dV:+.4f} e (should be ~0);  displaced charge 1/2 int|drho| = {0.5*np.abs(d).sum()*dV:.3f} e")
    print(f"drho range {d.min():.4f} .. {d.max():.4f} a.u.")

    # charge-displacement curve along centroid(B) -> centroid(A)
    xyz = ab["xyz"]
    cA, cB = xyz[:na].mean(0), xyz[na:].mean(0)
    u = (cA - cB) / np.linalg.norm(cA - cB)
    i, j, k = np.meshgrid(*[np.arange(m) for m in ab["n"]], indexing="ij")
    pts = ab["org"] + i[..., None] * ab["ax"][0] + j[..., None] * ab["ax"][1] + k[..., None] * ab["ax"][2]
    t = ((pts - cB) @ u) * BOHR  # Angstrom, 0 at centroid of B
    h = 0.05
    edges = np.arange(t.min(), t.max() + h, h)
    dq_slab, _ = np.histogram(t, edges, weights=d * dV)
    dq = np.cumsum(dq_slab)
    tc = edges[1:]
    # isodensity boundary: where rho_A = rho_B along the slabs
    rA, _ = np.histogram(t, edges, weights=a["data"])
    rB, _ = np.histogram(t, edges, weights=b["data"])
    L = np.linalg.norm(cA - cB) * BOHR
    win = (tc > 0) & (tc < L)
    tb = tc[win][np.argmin(np.abs(rA[win] - rB[win]))]
    ct = np.interp(tb, tc, dq)
    print(f"centroid distance {L:.2f} A; isodensity boundary at t = {tb:.2f} A; CT(boundary) = {ct:+.3f} e "
          f"({'A->B' if ct > 0 else 'B->A'})")
    np.savetxt(pre + "_cd.dat", np.c_[tc, dq], header="t/A (0 = centroid B)  Delta_q(t)/e", fmt="%.4f")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    ax.axhline(0, color="0.6", lw=0.7)
    ax.plot(tc, dq, color="k", lw=1.4)
    ax.axvline(tb, color="0.4", ls="--", lw=0.9)
    tA = (xyz[:na] - cB) @ u * BOHR
    tB = (xyz[na:] - cB) @ u * BOHR
    y0 = ax.get_ylim()[0]
    ax.plot(tB, np.full_like(tB, y0), "|", color="#21b0fe", ms=6, label=nameB)
    ax.plot(tA, np.full_like(tA, y0), "|", color="#fe218b", ms=6, label=nameA)
    ax.set_xlabel(f"position along {nameB}$\\rightarrow${nameA} axis (Å)")
    ax.set_ylabel("$\\Delta q(t)$ (e)")
    src, dst = (nameA, nameB) if ct > 0 else (nameB, nameA)
    ax.text(tb, ax.get_ylim()[1], f" {abs(ct):.3f} e  {src}$\\rightarrow${dst}", va="top", ha="left", fontsize=8)
    ax.text(0.02, 0.03, "$\\Delta q<0$: electrons moved\nleft $\\rightarrow$ right", transform=ax.transAxes,
            fontsize=6.5, color="0.35", va="bottom")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{pre}_cd.{ext}", dpi=300)


if __name__ == "__main__":
    main()
