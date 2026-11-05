import numpy as np
from sklearn.linear_model import LogisticRegression

from tte.weighting import stabilized_weights


def _confounded(n=6000, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3))
    a = rng.binomial(1, 1 / (1 + np.exp(-(0.2 + 0.9 * X[:, 0] - 0.5 * X[:, 1]))))
    ps = LogisticRegression().fit(X, a).predict_proba(X)[:, 1]
    return X, a, ps


def test_stabilized_weights_have_mean_about_one():
    X, a, ps = _confounded()
    w = stabilized_weights(ps, a, trunc=None)
    assert 0.95 < w.mean() < 1.05
    assert (w > 0).all()


def test_truncation_caps_extremes():
    X, a, ps = _confounded()
    raw = stabilized_weights(ps, a, trunc=None)
    cut = stabilized_weights(ps, a, trunc=(1, 99))
    assert cut.max() <= np.percentile(raw, 99) + 1e-12
    assert cut.min() >= np.percentile(raw, 1) - 1e-12


def test_weighting_removes_confounder_imbalance():
    X, a, ps = _confounded()
    w = stabilized_weights(ps, a, trunc=None)

    def smd(x, w):
        m1 = np.average(x[a == 1], weights=w[a == 1]); m0 = np.average(x[a == 0], weights=w[a == 0])
        return abs(m1 - m0) / x.std()

    assert smd(X[:, 0], np.ones(len(a))) > 0.3          # clearly imbalanced before weighting
    assert smd(X[:, 0], w) < 0.05                        # balanced after
