"""Retry a failed SLSQP solve with the same quadratic objective and constraints."""

import numpy as np
from types import SimpleNamespace
import rgmtm.error_budget as eb
import rgmtm.unified as unified
from contextlib import contextmanager
from . import conic_qp

original_frame = unified.solve_frame
original_minimize = eb.minimize
certificates = []


def conic_minimize(
    fun, x0, jac=None, method=None, bounds=None, constraints=(), options=None, **kwargs
):
    result = original_minimize(
        fun,
        x0,
        jac=jac,
        method=method,
        bounds=bounds,
        constraints=constraints,
        options=options,
        **kwargs,
    )
    if result.success or method != "SLSQP":
        return result
    x0 = np.asarray(x0, float)
    g0 = np.asarray(jac(x0), float)
    n = len(x0)
    step = max(1.0, float(np.max(abs(x0))))
    H = np.column_stack(
        [(np.asarray(jac(x0 + step * np.eye(n)[i])) - g0) / step for i in range(n)]
    )
    H = (H + H.T) / 2
    # Solve in increments about the supplied feasible initial point.
    A = []
    lo = []
    hi = []
    for c in ([constraints] if isinstance(constraints, dict) else constraints):
        M = np.asarray(c["jac"](x0), float)
        v = np.atleast_1d(c["fun"](x0))
        M = np.atleast_2d(M)
        A.extend(M)
        lo.extend(-v)
        hi.extend(-v if c["type"] == "eq" else np.full(len(v), np.inf))
    if bounds is not None:
        A.extend(np.eye(n))
        lo.extend(
            [(b[0] if b[0] is not None else -np.inf) - x for b, x in zip(bounds, x0)]
        )
        hi.extend(
            [(b[1] if b[1] is not None else np.inf) - x for b, x in zip(bounds, x0)]
        )
    delta, cert = conic_qp.solve(H, g0, np.asarray(A), np.asarray(lo), np.asarray(hi))
    x = x0 + delta
    # Verify that the Hessian reconstruction represents the same quadratic objective.
    expected = float(fun(x0) + g0 @ delta + 0.5 * delta @ H @ delta)
    actual = float(fun(x))
    assert abs(actual - expected) < 1e-7 * max(1, abs(actual)), (actual, expected)
    certificates.append(
        dict(
            original_error=str(result.message),
            dimension=n,
            certificate=cert,
            quadratic_value_residual=abs(actual - expected),
        )
    )
    return SimpleNamespace(
        success=True, x=x, fun=actual, message="same-QP Clarabel fallback, certified"
    )


def frame(*args, **kwargs):
    try:
        return original_frame(*args, **kwargs)
    except RuntimeError as ex:
        if "QP numerical failure" not in str(ex):
            raise
        print("RETRY SAME QP", str(ex), flush=True)
        eb.minimize = conic_minimize
        try:
            return original_frame(*args, **kwargs)
        finally:
            eb.minimize = original_minimize


@contextmanager
def retry_failed_frames():
    certificates.clear()
    unified.solve_frame = frame
    try:
        yield certificates
    finally:
        unified.solve_frame = original_frame
