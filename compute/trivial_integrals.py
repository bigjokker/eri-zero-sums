"""Trivial-zero integrals for Paper A, checked against their E_1 series.

Odd chi  (trivial zeros at the negative odd integers):
    int_y^inf dt / ((t^2 - 1) log t)   = sum_{k>=0} E_1((2k+1) log y)
Even chi (trivial zeros at the negative even integers):
    int_y^inf dt / (t(t^2 - 1) log t)  = sum_{k>=0} E_1((2k+2) log y)

Both follow from the geometric expansion of 1/(t^2-1) resp. 1/(t(t^2-1)) on t > 1
followed by u = log t.  The exponent is easy to get off by one, hence this check.
"""

import mpmath as mp

mp.mp.dps = 30
YS = [2, 3, 10, 100]


def check(y):
    y = mp.mpf(y)
    L = mp.log(y)
    odd_int = mp.quad(lambda t: 1 / ((t ** 2 - 1) * mp.log(t)), [y, mp.inf])
    odd_ser = mp.nsum(lambda k: mp.e1((2 * k + 1) * L), [0, mp.inf])
    even_int = mp.quad(lambda t: 1 / (t * (t ** 2 - 1) * mp.log(t)), [y, mp.inf])
    even_ser = mp.nsum(lambda k: mp.e1((2 * k + 2) * L), [0, mp.inf])
    return odd_int, odd_ser, even_int, even_ser


if __name__ == "__main__":
    worst = mp.mpf(0)
    lines = [f"mp.dps = {mp.mp.dps}", ""]
    for y in YS:
        oi, os_, ei, es = check(y)
        do, de = abs(oi - os_), abs(ei - es)
        worst = max(worst, do, de)
        lines.append(f"y = {y}")
        lines.append(f"  odd   integral {mp.nstr(oi, 20)}")
        lines.append(f"        series   {mp.nstr(os_, 20)}   diff {mp.nstr(do, 3)}")
        lines.append(f"  even  integral {mp.nstr(ei, 20)}")
        lines.append(f"        series   {mp.nstr(es, 20)}   diff {mp.nstr(de, 3)}")
    lines.append("")
    lines.append(f"worst discrepancy {mp.nstr(worst, 3)}")
    out = "\n".join(lines) + "\n"
    with open("data/trivial_integrals.txt", "w", encoding="utf-8") as f:
        f.write(out)
    print(out)
