"""Stabilized inverse-probability-of-treatment weights."""
import numpy as np


def stabilized_weights(ps, treated, trunc=(1, 99)):
    """SW_i = P(A = a_i) / P(A = a_i | L_i), optionally truncated at the given percentiles."""
    ps = np.asarray(ps, float)
    a = np.asarray(treated)
    p = a.mean()
    w = np.where(a == 1, p / ps, (1 - p) / (1 - ps))
    if trunc is not None:
        lo, hi = np.percentile(w, trunc)
        w = np.clip(w, lo, hi)
    return w
