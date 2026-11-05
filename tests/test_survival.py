import numpy as np
from lifelines import KaplanMeierFitter

from tte.survival import weighted_aj, at


def _data(n=800, seed=0, competing=True):
    rng = np.random.default_rng(seed)
    t1 = rng.exponential(60, n)                                   # time to CVD death
    t2 = rng.exponential(90, n) if competing else np.full(n, np.inf)   # time to other death
    c = rng.uniform(5, 150, n)                                    # censoring
    first = np.minimum(t1, t2)
    ev = np.where(first <= c, np.where(t1 < t2, 1, 2), 0)
    t = np.ceil(np.minimum(first, c))                             # integer months -> ties, as in NHANES
    return t, ev


def test_cif_equals_one_minus_km_without_competing_events():
    t, ev = _data(competing=False)
    u, S, c1, c2 = weighted_aj(t, ev, np.ones(len(t)))
    km = KaplanMeierFitter().fit(t, ev == 1)
    for x in (12, 36, 60, 120):
        assert abs(at(u, c1, x) - (1 - km.predict(x))) < 1e-9
    assert np.allclose(c2, 0)


def test_probabilities_sum_to_one_with_weights():
    t, ev = _data()
    w = np.random.default_rng(1).uniform(0.3, 3.0, len(t))
    u, S, c1, c2 = weighted_aj(t, ev, w)
    assert np.allclose(S + c1 + c2, 1.0)
    assert np.all(np.diff(c1) >= -1e-12) and np.all(np.diff(c2) >= -1e-12)   # CIFs never decrease


def test_integer_weights_equal_row_duplication():
    t, ev = _data(n=300)
    w = np.random.default_rng(2).integers(1, 4, len(t))
    a = weighted_aj(t, ev, w)
    b = weighted_aj(np.repeat(t, w), np.repeat(ev, w), np.ones(w.sum()))
    for x, y in zip(a, b):
        assert np.allclose(x, y)
