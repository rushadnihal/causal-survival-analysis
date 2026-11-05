"""Weighted Aalen-Johansen cumulative incidence with one competing event."""
import numpy as np


def weighted_aj(t, ev, w):
    """Weighted Aalen-Johansen estimator.

    t  : follow-up times
    ev : 0 = censored, 1 = event of interest (CVD death), 2 = competing event (other death)
    w  : non-negative weights (e.g. stabilized IPTW, optionally x survey weight)

    Returns (times, S, CIF1, CIF2): step functions evaluated at the distinct observed times.
    Tied events are handled by aggregating weights at each distinct time; at each time the
    risk set contains everyone with t >= that time (events precede censoring at ties).
    """
    t, ev, w = np.asarray(t, float), np.asarray(ev), np.asarray(w, float)
    o = np.argsort(t, kind="stable")
    t, ev, w = t[o], ev[o], w[o]
    u, idx = np.unique(t, return_index=True)
    at_risk = np.cumsum(w[::-1])[::-1][idx]
    d1 = np.add.reduceat(w * (ev == 1), idx)
    d2 = np.add.reduceat(w * (ev == 2), idx)
    S, c1, c2 = 1.0, 0.0, 0.0
    S_out, c1_out, c2_out = [], [], []
    for n, a, b in zip(at_risk, d1, d2):
        c1 += S * a / n
        c2 += S * b / n
        S *= 1 - (a + b) / n
        S_out.append(S)
        c1_out.append(c1)
        c2_out.append(c2)
    return u, np.array(S_out), np.array(c1_out), np.array(c2_out)


def at(times, vals, x):
    """Value of a right-continuous step function (times, vals) at x; 0 before the first time."""
    i = np.searchsorted(times, x, side="right") - 1
    return 0.0 if i < 0 else vals[i]


def restricted_mean(times, vals, tau, initial):
    """Area under a right-continuous step function on [0, tau].

    `initial` is the function's value before the first jump (1.0 for a survival curve,
    0.0 for a cumulative incidence curve). For S this is the restricted mean survival time;
    for CIF_k it is the mean time lost to cause k before tau.
    """
    times, vals = np.asarray(times, float), np.asarray(vals, float)
    keep = times < tau
    pts = np.concatenate(([0.0], times[keep], [tau]))
    heights = np.concatenate(([initial], vals[keep]))
    return float(np.sum(heights * np.diff(pts)))
