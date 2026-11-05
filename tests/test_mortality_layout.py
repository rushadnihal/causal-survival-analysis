"""Regression tests for the NCHS fixed-width mortality layout (bugs found during development:
PERMTH columns were off by one, and 6-digit 2017-18 SEQNs were truncated to 5 characters)."""
from tte.mortality import read_mortality


def _line(seqn, elig, mort, ucod, diab, hyp, pint, pexm):
    ch = [" "] * 48

    def put(s, start_col):                  # 1-based start column, left-justified
        for i, c in enumerate(str(s)):
            ch[start_col - 1 + i] = c

    put(seqn, 1); put(elig, 15); put(mort, 16)
    put(ucod, 17); put(diab, 20); put(hyp, 21)
    put(pint, 43); put(pexm, 46)
    return "".join(ch).rstrip() + "\n"


def test_fields_parse_at_documented_columns(tmp_path):
    p = tmp_path / "mort.dat"
    p.write_text(_line(31130, 1, 1, "004", 0, 1, 10, 8) +        # CVD-irrelevant death (accident), short follow-up
                 _line(93703, 1, 0, "", "", "", 159, 158))        # alive, 6-digit SEQN, long follow-up
    m = read_mortality(str(p))
    assert list(m.SEQN) == [31130, 93703]
    assert list(m.MORTSTAT) == [1, 0]
    assert m.UCOD_LEADING.iloc[0] == 4 and m.UCOD_LEADING.isna().iloc[1]
    assert list(m.PERMTH_INT) == [10, 159]
    assert list(m.PERMTH_EXM) == [8, 158]
