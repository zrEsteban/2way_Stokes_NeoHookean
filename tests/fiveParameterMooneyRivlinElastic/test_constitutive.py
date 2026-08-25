#!/usr/bin/env python3
"""Dependency-free constitutive verification for the compiled solids4foam law.

The equations here deliberately mirror fiveParameterMooneyRivlinElastic.C and
write all pointwise results to CSV for auditability.
"""

import csv
import math
from pathlib import Path

LAMBDAS = (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4)
MATERIALS = {
    "young": (1.7668e6, -1.3777e6, 5.0640e5, -4.6800e4, 4.1000e3),
    "aged": (1.6240e5, -7.2200e4, 3.0740e5, -1.1830e5, 1.3400e4),
}


def invariants(lam):
    return lam * lam + 2.0 / lam, 2.0 * lam + 1.0 / (lam * lam)


def energy(lam, c):
    i1, i2 = invariants(lam)
    x, y = i1 - 3.0, i2 - 3.0
    c10, c01, c20, c11, c02 = c
    return c10*x + c01*y + c20*x*x + c11*x*y + c02*y*y


def stress_difference(lam, c):
    i1, i2 = invariants(lam)
    x, y = i1 - 3.0, i2 - 3.0
    c10, c01, c20, c11, c02 = c
    w1 = c10 + 2*c20*x + c11*y
    w2 = c01 + c11*x + 2*c02*y
    # sigma_11 - sigma_22; hydrostatic pressure and dev() cancel in difference.
    return 2*w1*(lam*lam - 1/lam) - 2*w2*(1/(lam*lam) - lam)


def central_derivative(fun, x, h):
    return (fun(x + h) - fun(x - h))/(2*h)


def analytic_energy_derivative(lam, c):
    i1, i2 = invariants(lam)
    x, y = i1 - 3.0, i2 - 3.0
    c10, c01, c20, c11, c02 = c
    w1 = c10 + 2*c20*x + c11*y
    w2 = c01 + c11*x + 2*c02*y
    di1 = 2*lam - 2/(lam*lam)
    di2 = 2 - 2/(lam**3)
    return w1*di1 + w2*di2


def relative_error(a, b, floor=1.0):
    return abs(a-b)/max(abs(a), abs(b), floor)


def main():
    out_dir = Path(__file__).resolve().parent / "results"
    out_dir.mkdir(exist_ok=True)
    failures = []
    rows = []

    for name, coeffs in MATERIALS.items():
        mu0 = 2*(coeffs[0] + coeffs[1])
        expected = 7.782e5 if name == "young" else 1.804e5
        if relative_error(mu0, expected) > 1e-12:
            failures.append(f"{name}: mu0={mu0}, expected={expected}")

        if abs(energy(1.0, coeffs)) > 1e-12:
            failures.append(f"{name}: undeformed energy is not zero")
        if abs(stress_difference(1.0, coeffs)) > 1e-9:
            failures.append(f"{name}: undeformed stress difference is not zero")

        for lam in LAMBDAS:
            h = 1e-5*max(1.0, abs(lam))
            w = energy(lam, coeffs)
            sigma = stress_difference(lam, coeffs)
            d_w_analytic = analytic_energy_derivative(lam, coeffs)
            d_w_fd = central_derivative(lambda q: energy(q, coeffs), lam, h)
            sigma_from_w = lam*d_w_fd
            tangent_analytic = central_derivative(
                lambda q: q*analytic_energy_derivative(q, coeffs), lam, h/10
            )
            tangent_fd = central_derivative(
                lambda q: stress_difference(q, coeffs), lam, h
            )
            e_energy = relative_error(d_w_analytic, d_w_fd)
            e_stress = relative_error(sigma, sigma_from_w)
            e_tangent = relative_error(tangent_analytic, tangent_fd)
            rows.append((name, lam, w, sigma, sigma_from_w, tangent_analytic,
                         tangent_fd, e_energy, e_stress, e_tangent))
            derivative_scale = max(abs(d_w_analytic), abs(d_w_fd))
            stress_scale = max(abs(sigma), abs(sigma_from_w))
            energy_failed = derivative_scale > 100 and e_energy > 2e-8
            stress_failed = stress_scale > 100 and e_stress > 2e-8
            tangent_failed = e_tangent > 1e-7
            if energy_failed or stress_failed or tangent_failed:
                failures.append(
                    f"{name} lambda={lam}: errors "
                    f"energy={e_energy:g}, stress={e_stress:g}, tangent={e_tangent:g}"
                )

    # Exact algebraic regression against the installed three-parameter law
    # with c11=0 (the conventional two-parameter specialization).
    for lam in LAMBDAS:
        c = (2.3e5, 1.1e5, 0.0, 0.0, 0.0)
        i1, i2 = invariants(lam)
        old = 2*c[0]*(lam*lam - 1/lam) - 2*c[1]*(1/(lam*lam) - lam)
        new = stress_difference(lam, c)
        if old != new:
            failures.append(f"two-parameter regression failed at lambda={lam}")

    eps = 1e-6
    young_tangent = central_derivative(
        lambda q: stress_difference(q, MATERIALS["young"]), 1.0, eps
    )
    aged_tangent = central_derivative(
        lambda q: stress_difference(q, MATERIALS["aged"]), 1.0, eps
    )
    if not young_tangent > aged_tangent > 0:
        failures.append("young PDMS is not initially stiffer than aged PDMS")

    with (out_dir / "uniaxial.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("material", "lambda", "W_Pa", "sigmaDifference_Pa",
                         "sigmaFromEnergy_Pa", "analyticTangent_Pa",
                         "finiteDifferenceTangent_Pa", "energyRelativeError",
                         "stressRelativeError", "tangentRelativeError"))
        writer.writerows(rows)

    max_errors = [max(row[i] for row in rows) for i in (7, 8, 9)]
    report = (
        f"tests={'FAIL' if failures else 'PASS'}\n"
        f"youngMu0={2*(MATERIALS['young'][0]+MATERIALS['young'][1]):.12g}\n"
        f"agedMu0={2*(MATERIALS['aged'][0]+MATERIALS['aged'][1]):.12g}\n"
        f"youngInitialTangent={young_tangent:.12g}\n"
        f"agedInitialTangent={aged_tangent:.12g}\n"
        f"maxEnergyRelativeError={max_errors[0]:.12g}\n"
        f"maxStressRelativeError={max_errors[1]:.12g}\n"
        f"maxTangentRelativeError={max_errors[2]:.12g}\n"
    )
    (out_dir / "summary.txt").write_text(report)
    print(report, end="")
    if failures:
        print("\n".join(failures))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
