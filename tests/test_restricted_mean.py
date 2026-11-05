import numpy as np

from tte.survival import weighted_aj, restricted_mean


def _data(n=800, seed=0):
    rng = np.random.default_rng(seed)
    t1 = rng.exponential(60, n)            # time to CVD death
    t2 = rng.exponential(90, n)            # time to other death
    c = rng.uniform(5, 150, n)             # censoring
    first = np.minimum(t1, t2)
    ev = np.where(first <= c, np.where(t1 < t2, 1, 2), 0)
    t = np.ceil(np.minimum(first, c))
    return t, ev


def test_restricted_mean_hand_computed():
    times, S, F = [1.0, 2.0], [0.5, 0.25], [0.5, 0.75]
    assert np.isclose(restricted_mean(times, S, 3.0, 1.0), 1 * 1 + 0.5 * 1 + 0.25 * 1)   # 1.75
    assert np.isclose(restricted_mean(times, F, 3.0, 0.0), 0 * 1 + 0.5 * 1 + 0.75 * 1)   # 1.25
    assert np.isclose(restricted_mean(times, S, 1.5, 1.0), 1.0 + 0.5 * 0.5)               # horizon inside an interval


def test_mean_times_add_up_to_horizon():
    t, ev = _data()
    w = np.random.default_rng(3).uniform(0.5, 2.0, len(t))
    u, S, c1, c2 = weighted_aj(t, ev, w)
    tau = 60.0
    total = restricted_mean(u, S, tau, 1.0) + restricted_mean(u, c1, tau, 0.0) + restricted_mean(u, c2, tau, 0.0)
    assert np.isclose(total, tau)
