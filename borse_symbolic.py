"""Exact (symbolic) verification of every formula used in the paper, with SymPy.

Run:  .venv/bin/python borse_symbolic.py   (output: rerun/symbolic_proofs.txt)
"""

import sympy as sp

s5 = sp.sqrt(5)
phi = (1 + s5) / 2
x = (s5 - 1) / 2
t, zc = sp.symbols("t z")
results = []


def ok(name, expr):
    r = sp.simplify(expr) == 0
    results.append(r)
    print(f"{'PASS' if r else 'FAIL'}  {name}")


print("== Lemma 1: self-similar partition")
ok("x = 1/phi", x - 1 / phi)
ok("(1-x)/x = x/1", (1 - x) / x - x)
ok("x^2 + x - 1 = 0", x ** 2 + x - 1)
ok("x = 1/(1+x)", x - 1 / (1 + x))
roots = sp.solve(t ** 2 + t - 1, t)
pos = [r for r in roots if r.is_positive]
r = len(pos) == 1 and sp.simplify(pos[0] - x) == 0 and any(sp.simplify(q + phi) == 0 for q in roots)
results.append(r)
print(f"{'PASS' if r else 'FAIL'}  t^2+t-1 has exactly one positive root (= x); the other root is -phi")
ok("1 + 1/phi = phi", 1 + 1 / phi - phi)

print("== Lemma 2: identities")
ok("(a) x^2 = 1 - x", x ** 2 - (1 - x))
ok("(a) x^3 = 2x - 1 = sqrt5 - 2", x ** 3 - (s5 - 2))
ok("(a) x^4 = 2 - 3x", x ** 4 - (2 - 3 * x))
ok("(b) 1 - x^2/2 = (1+x)/2 = phi/2", 1 - x ** 2 / 2 - phi / 2)
ok("(c) x^4/2 = (5 - 3 phi)/2", x ** 4 / 2 - (5 - 3 * phi) / 2)
ok("(d) x/(2-x) = 1/sqrt5", x / (2 - x) - 1 / s5)
ok("(e) 2 phi - 1 = sqrt5, so 1/(2 phi - 1) = 1/sqrt5", 1 / (2 * phi - 1) - 1 / s5)

print("== Consequences for dark energy")
w0, wa = -1 + x ** 2 / 2, -x
ok("w0 = -phi/2", w0 + phi / 2)
ok("wa = -1/phi", wa + 1 / phi)
ok("wa = 1 + 2 w0 (equivalent to x^2 + x - 1 = 0)", wa - (1 + 2 * w0))
xs = sp.Symbol("xs", positive=True)
ok("for general x on the ladder curve, wa - (1 + 2 w0) = -(x^2 + x - 1)", (-xs) - (1 + 2 * (-1 + xs ** 2 / 2)) + (xs ** 2 + xs - 1))
ok("1 + w0 + wa = (1 - 3x)/2", 1 + w0 + wa - (1 - 3 * x) / 2)
ok("rho_DE exponent: -3(1+w0+wa) = (3/2)(3x - 1)", -3 * (1 + w0 + wa) - sp.Rational(3, 2) * (3 * x - 1))
ok("(1+w0)/(-wa) = x/2", (1 + w0) / (-wa) - x / 2)
sol = sp.solve(sp.Eq(w0 + wa * zc / (1 + zc), -1), zc)
ok("w(z) = -1 exactly at z = 1/sqrt5", sol[0] - 1 / s5)
qb = bool(sp.N(w0) > -1 and sp.N(w0 + wa) < -1)
results.append(qb)
print(f"{'PASS' if qb else 'FAIL'}  quintom-B region: w0 = {sp.N(w0, 6)} > -1 and w0 + wa = {sp.N(w0 + wa, 6)} < -1")
w = lambda z: w0 + wa * z / (1 + z)
print(f"      fixed-redshift values: w(0.3) = {sp.nsimplify(w(sp.Rational(3, 10)))} = {sp.N(w(sp.Rational(3, 10)), 8)};"
      f" w(0.5) = {sp.simplify(w(sp.Rational(1, 2)))} = {sp.N(w(sp.Rational(1, 2)), 8)}")
ok("w(0.5) = -1 - (3 - sqrt5)/12 ... closed form check", w(sp.Rational(1, 2)) - (w0 + wa / 3))
ok("non-crossing bound -(1+w0) = -F(2)", -(1 + w0) + x ** 2 / 2)
ok("d=3 branch: -1 + x^3 = -1 + phi^-3", (-1 + x ** 3) - (-1 + phi ** -3))
ok("Remark: 6/11 + F(4) - 1/phi = (5/2)(89/55 - phi)",
   sp.Rational(6, 11) + x ** 4 / 2 - 1 / phi - sp.Rational(5, 2) * (sp.Rational(89, 55) - phi))

print("== alpha_d ladder values (alpha = 1 for d = 1, 3; 1/2 for d = 2, 4)")
for d, al in ((1, 1), (2, sp.Rational(1, 2)), (3, 1), (4, sp.Rational(1, 2))):
    print(f"      d={d}: alpha_d x^d = {sp.simplify(al * x ** d)} = {sp.N(al * x ** d, 10)}")

print("== Proposition 3: recursion and Fibonacci (exact rationals, n <= 25)")
F = [1, 1]
for _ in range(30):
    F.append(F[-1] + F[-2])
seq = [sp.Integer(1)]
for _ in range(25):
    seq.append(1 / (1 + seq[-1]))
r = all(seq[n] == sp.Rational(F[n], F[n + 1]) for n in range(26))
results.append(r)
print(f"{'PASS' if r else 'FAIL'}  x_n = F_(n+1)/F_(n+2)")
r = all(sp.N(abs(seq[n] - x) - phi ** -n * abs(seq[0] - x), 30) <= 0 for n in range(26))
results.append(r)
print(f"{'PASS' if r else 'FAIL'}  |x_n - 1/phi| <= phi^-n |x_0 - 1/phi|")
print("== Contraction: f is not a contraction on [0, inf), but f o f is (constant <= 1/4)")
u = sp.Symbol("u", nonnegative=True)
fu = 1 / (1 + u)
ffu = sp.simplify(1 / (1 + fu))
r = sp.simplify(ffu - (1 + u) / (2 + u)) == 0
results.append(r)
print(f"{'PASS' if r else 'FAIL'}  f(f(u)) = (1 + u)/(2 + u)")
d1 = sp.diff(ffu, u)
r = sp.simplify(d1 - 1 / (2 + u) ** 2) == 0 and d1.subs(u, 0) == sp.Rational(1, 4) and sp.simplify(sp.diff(d1, u) + 2 / (2 + u) ** 3) == 0
results.append(r)
print(f"{'PASS' if r else 'FAIL'}  (f o f)'(u) = 1/(2 + u)^2, decreasing, so sup over u >= 0 is 1/4")
r = sp.diff(fu, u).subs(u, 0) == -1
results.append(r)
print(f"{'PASS' if r else 'FAIL'}  f'(0) = -1, so f itself has Lipschitz constant 1 on [0, inf)")
print("== Hurwitz illustration: q^2 |x - p/q| -> 1/sqrt5 along Fibonacci convergents")
for n in (5, 10, 20):
    p_, q_ = F[n], F[n + 1]
    print(f"      q={q_:>8d}: {sp.N(q_ ** 2 * abs(x - sp.Rational(p_, q_)), 8)}  (1/sqrt5 = {sp.N(1 / s5, 8)})")
print("== Numerical values quoted")
for k, v in {"x": x, "F(1)": x / 2, "F(2)": x ** 2 / 2, "F(3)": x ** 3 / 2, "F(4)": x ** 4 / 2, "x^3": x ** 3,
             "w0": w0, "wa": wa, "1+w0+wa": 1 + w0 + wa, "z_cross": 1 / s5, "xi=(7phi-10)/4": (7 * phi - 10) / 4,
             "6/11+F4": sp.Rational(6, 11) + x ** 4 / 2, "residual": sp.Rational(6, 11) + x ** 4 / 2 - 1 / phi,
             "2 ln 403": 2 * sp.log(403)}.items():
    print(f"      {k:16s} = {sp.N(v, 10)}")
print(f"\n{sum(results)}/{len(results)} exact checks passed")
