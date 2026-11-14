#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the throat and the potentials of the noncommutative wormhole.

    f(x) = 1 - erf(mu x)/x + 2 mu exp(-mu^2 x^2)/sqrt(pi)       (as in matrixassembler_bc1.mpl)
    V(z) = l(l + 1)/z^2 + (1 - s)(1 + 2s)/(2z) df/dz,   z = x/x0 > 1

The throat x0 is the larger zero of f, found as the assemblers find it; the extremal mass solves
f = 0 and df/dx = 0 together. From this folder,
    python3 plot_readme_figure.py
writes ncwh-throat.png (Python 3 with NumPy, SciPy and Matplotlib). The style is the stamped copy
manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import brentq, fsolve
from scipy.special import erf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402


def f(x, mu):
    return 1 - erf(mu * x) / x + 2 * mu * np.exp(-(mu * x) ** 2) / np.sqrt(np.pi)


def dfdx(x, mu):
    return erf(mu * x) / x**2 - 2 * mu * np.exp(-(mu * x) ** 2) / (np.sqrt(np.pi) * x) \
        - 4 * mu**3 * x * np.exp(-(mu * x) ** 2) / np.sqrt(np.pi)


def throat(mu):
    """The larger zero x0 = x+ of f, beyond the minimum of f."""
    xmin = brentq(lambda t: dfdx(t, mu), 0.3, 1.5)
    return brentq(lambda t: f(t, mu), xmin, 3.0, xtol=1e-14)


xe, mu_e = fsolve(lambda v: [f(v[0], v[1]), dfdx(v[0], v[1])], (0.79, 1.90), xtol=1e-14)
print(f"extremal configuration: mu_e = {mu_e:.10f}, x_e = {xe:.10f}")

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.75))

# (a) the function f and its larger zero, the throat x0.
x = np.linspace(1e-3, 4.0, 4000)
for k, (mu, words) in enumerate([(mu_e, "extremal"), (1.91, "near-extremal"), (2.3, "")]):
    style = hs.series(k)
    style["marker"] = None
    label = rf"$\mu = {mu:.6g}$" + (f" ({words})" if words else "")
    left.plot(x, f(x, mu), lw=1.1, label=label, **style)
    x0 = xe if k == 0 else throat(mu)
    left.plot([x0], [0.0], linestyle="none", marker="o", markersize=4.2,
              markerfacecolor="white", markeredgecolor=hs.SERIES[k], markeredgewidth=0.9)
    print(f"mu = {mu:.10g}: throat x0 = {x0:.10f}")
left.plot(x[x > 0.5], 1 - 1 / x[x > 0.5], color=hs.SLATE, lw=0.8, linestyle=(0, (5.0, 1.5, 1.5, 1.5)),
          label=r"$1 - 1/x$ (classical)")
left.axhline(0.0, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
left.set_xlim(0, 4.0)
left.set_ylim(-0.3, 1.05)
left.set_xlabel(r"$x = r/(2M)$")
left.set_ylabel(r"$f(x)$")
hs.panel(left, r"(a) the function $f$; its larger zero is the throat")
hs.legend_outside(left, where="below", ncol=2)

# (b) the potentials of the three perturbations, l = 2, near extremality (mu = 1.91).
mu = 1.91
x0 = throat(mu)
z = np.linspace(1.0, 6.0, 2000)
ell = 2
for k, s in enumerate((0, 1, 2)):
    style = hs.series(k)
    style["marker"] = None
    V = ell * (ell + 1) / z**2 + (1 - s) * (1 + 2 * s) / (2 * z) * x0 * dfdx(x0 * z, mu)
    right.plot(z, V, lw=1.1, label=rf"$s = {s}$", **style)
    print(f"mu = {mu}, l = {ell}, s = {s}: V(1) = {V[0]:.6f}, max V = {V.max():.6f} at z = {z[V.argmax()]:.4f}")
right.set_xlim(1.0, 6.0)
right.set_ylim(0.0, 6.5)
right.set_xlabel(r"$z = x/x_0$")
right.set_ylabel(r"$V_{s,\ell}(z)$")
hs.panel(right, r"(b) potentials for $\ell = 2$ at $\mu = 1.91$")
hs.legend_outside(right, where="below", ncol=3)

out = HERE / "ncwh-throat.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
