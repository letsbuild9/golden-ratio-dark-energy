"""Automated check that every computed number quoted in paper/paper.tex matches the analysis outputs (version 2).

Reads the outputs in rerun/ (produced by a clean run of all scripts), formats each value exactly as the paper
prints it, and checks that the formatted string occurs in the LaTeX source. Prints PASS/FAIL per item.

Run:  .venv/bin/python verify_paper.py
"""

import re
import sys

TEX = open("paper/paper.tex").read()
R = lambda f: open(f"rerun/{f}").read()
results = []


def check(label, needle):
    ok = needle in TEX
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {label:55s} -> {needle}")


def rng(vals, fmt="{:.1f}", sep="$--$"):
    return fmt.format(min(vals)) + sep + fmt.format(max(vals))


f1 = lambda v: f"{v:.1f}"

# ---------------------------------------------------------------- Table 1 (T6)
v2 = R("v2_results.txt")
t6 = v2.split("T3.")[0]
rows = re.findall(r"^(DESI\+CMB\S*)\s+(LCDM|Borse d=2|Borse d=3|w0waCDM)\s+\d\s+([\d.]+)\s+\d+\s+([\d.]+)\s+[\d.]+\s+([\d.]+)\s+([\d.]+)\s+(\S+)$", t6, re.M)
lab = {"DESI+CMB": "DESI+CMB", "DESI+CMB+PantheonPlus": "+Pantheon+", "DESI+CMB+Union3": "+Union3",
       "DESI+CMB+DESY5": "+DESY5", "DESI+CMB+DESDovekie": "+DES-Dovekie"}
tex_model = {"LCDM": r"\LCDM", "Borse d=2": r"$\gr$-ladder", "w0waCDM": r"$w_0w_a$"}
T = {}
for data, model, c2, red, aic, bic, ns in rows:
    T[(data, model)] = dict(c2=float(c2), aic=aic, bic=bic, ns=ns, red=red)
    if model in tex_model:
        a, b = (rf"\textbf{{{aic}}}", rf"\textbf{{{bic}}}") if model == "Borse d=2" else (aic, bic)
        prefix = f"{lab[data]} & " if model == "LCDM" else " & "
        check(f"Table 1 {lab[data]} {model}", f"{prefix}{tex_model[model]} & {c2} & {red} & {a} & {b} & {ns if ns != '--' else '--'}")
datas = list(lab)
ladder = [float(T[(d, 'Borse d=2')]['ns']) for d in datas]
lcdm = [float(T[(d, 'LCDM')]['ns']) for d in datas]
d3 = [float(T[(d, 'Borse d=3')]['ns']) for d in datas]
check("phi-ladder vs free range (abstract and text)", f"lies ${rng(ladder)}\\sigma$ from the free fit")
check("phi-ladder vs free range (conclusion)", f"lies ${rng(ladder)}\\sigma$ from the best free fit")
check("LCDM disfavored range", f"\\LCDM\\ is disfavored at ${rng(lcdm)}\\sigma$")
check("d=3 branch range (Table 1 caption)", f"The $d = 3$ branch is ${rng(d3)}\\sigma$")
dchi = [T[(d, 'LCDM')]['c2'] - T[(d, 'Borse d=2')]['c2'] for d in datas]
check("phi-ladder preferred over LCDM by dchi2", f"$\\Delta\\chi^2 = {rng(dchi)}$")
dov = re.search(r"DESI\+CMB\+DESDovekie.*\n.*free best fit: w0 = (-[\d.]+), wa = (-[\d.]+)", t6)
check("DES-Dovekie free best fit", f"$(w_0, w_a) = ({dov.group(1)}, {dov.group(2)})$")
fit_txt = R("fit_results.txt")
h0 = [float(x) for x in re.findall(r"^Borse d=2\s+3\s+[\d.]+\s+-?[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s+[\d.]+", fit_txt, re.M)]
om = [float(x) for x in re.findall(r"^Borse d=2\s+3\s+[\d.]+\s+-?[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)", fit_txt, re.M)]
check("phi-ladder backgrounds H0 and Om", f"$H_0 \\simeq {rng(h0)}$\\,km\\,s$^{{-1}}$\\,Mpc$^{{-1}}$ and $\\Om \\simeq {rng(om, '{:.3f}')}$")

# ---------------------------------------------------------------- look-elsewhere (T3)
t3 = v2.split("T3.")[1].split("T1/T2.")[0]
ranks = [int(x) for x in re.findall(r"rank\s+(\d+)/403", t3)]
check("look-elsewhere ranks", f"ranks {ranks[0]}th (Pantheon+), {ranks[1]}th (Union3), {ranks[2]}th (DESY5) and {ranks[3]}nd (DES-Dovekie) out of 403")
cnt = [int(x) for x in re.findall(r"within 1 sigma of free fit:\s+(\d+)", t3)]
check("pairs within 1 sigma per sample", f"{min(cnt)}--{max(cnt)} pairs from the family")
only = re.findall(r"within 1 sigma for all of \[[^\]]+\]: (\[.*\])", t3)
ok = all(o == "[(-0.809, -0.618)]" for o in only) and len(only) == 2
results.append(ok)
print(f"{'PASS' if ok else 'FAIL'}  only the phi pair within 1 sigma for all three samples (both DES versions) -> {only}")

# ---------------------------------------------------------------- free x (T1/T2)
t1 = v2.split("T1/T2.")[1].split("T1b/T8")[0]
fx = re.findall(r"^(DESI\+CMB\S*): best x = ([\d.]+)\s+\(68% interval on 0.01 grid: \[([\d.]+), ([\d.]+)\]\);\s+1/phi is at dchi2 = ([\d.]+)", t1, re.M)
for d, xb, lo, hi, dphi in fx:
    check(f"free-x {lab[d]}", f"{lab[d]} & {xb} & [{lo}, {hi}]")
xbs = [float(r[1]) for r in fx]
check("free-x range (abstract and text)", f"$x = {rng(xbs, '{:.2f}')}$")
los, his = max(float(r[2]) for r in fx), min(float(r[3]) for r in fx)
check("overlap of the five 68% ranges", f"overlap only in $[{los:.2f}, {his:.2f}]$")
pp = [r for r in fx if r[0] == "DESI+CMB+PantheonPlus"][0]
check("1/phi at the edge for Pantheon+", f"($\\Delta\\chi^2 = {pp[4]}$)")
named = {}
for d, line in re.findall(r"^(DESI\+CMB\S*): best x.*\n\s+dchi2 at named constants: (.*)$", t1, re.M):
    named[d] = dict(kv.split("=") for kv in line.split(", "))
inside_all = [k for k in named["DESI+CMB"] if all(float(named[d][k]) < 1 for d in named)]
ok = inside_all == ["1/phi"]
results.append(ok)
print(f"{'PASS' if ok else 'FAIL'}  1/phi is the only named constant inside every 68% range -> {inside_all}")

# ---------------------------------------------------------------- official chains (T1b/T8/T4/T5)
tb = v2.split("T1b/T8/T4/T5.")[1].split("T12.")[0]
off = {}
for line in tb.splitlines():
    m = re.match(r"^(.{32}) +([\d.]+)%\s+([\d.]+) \| +([\d.]+) \[([\d.]+),([\d.]+)\] \|\s+(\S+)\s+(\S+).*?\|\s+([\d.]+)\s+(-[\d.]+)\+-([\d.]+)\s+(-[\d.]+)\s+(-?[\d.]+) \|\s+(-[\d.]+)\+-([\d.]+)\s+(-?[\d.]+) \|\s+(-[\d.]+)\+-([\d.]+)\s+(-?[\d.]+)", line)
    if m:
        off[m.group(1).strip()] = m.groups()
names = {"DESI+CMB (CamSpec) [baseline]": "DESI+CMB (baseline)", "+Pantheon+ [baseline]": "+Pantheon+", "+Union3 [baseline]": "+Union3",
         "+DESY5 [baseline]": "+DESY5", "DESI+CMB (Planck PR3 plik)": "DESI+CMB, Planck PR3 plik", "+DESY5 (plik)": "+DESY5, Planck PR3 plik",
         "+DESY5 (HiLLiPoP/LoLLiPoP)": "+DESY5, HiLLiPoP/LoLLiPoP", "+DESY5 with z>0.1": "+DESY5 with $z > 0.1$",
         "DESI BAO only (no CMB)": "DESI BAO only (no CMB)", "DESI BAO + DESY5 (no CMB)": "DESI BAO + DESY5 (no CMB)"}
for k, tex in names.items():
    g = off[k]
    check(f"Table 2 {tex}", f"{tex} & {float(g[1]):.1f}\\% & {g[2]}")
alt = [float(off[k][2]) for k in ("DESI+CMB (Planck PR3 plik)", "+DESY5 (plik)", "+DESY5 (HiLLiPoP/LoLLiPoP)")]
check("alternative CMB likelihood range", f"(${rng(alt, '{:.2f}')}\\sigma$)")
nocmb = [float(off[k][2]) for k in ("DESI BAO only (no CMB)", "DESI BAO + DESY5 (no CMB)")]
check("no-CMB range", f"(${rng(nocmb, '{:.2f}')}\\sigma$)")
check("low-z SN cut", f"low-redshift SN cut (${float(off['+DESY5 with z>0.1'][2]):.2f}\\sigma$)")
curve = [float(off[k][3]) for k in list(names)[:4]]
check("official-chain curve x range", f"(${rng(curve, '{:.3f}')}$)")
for k, tex in list(names.items())[:4]:
    g = off[k]
    zp, wp, sw, wt, dev = g[8], g[9], g[10], g[11], float(g[12]) + 0.0  # -0.00 prints as 0.00 in the paper
    tex0 = "DESI+CMB" if k.startswith("DESI+CMB") else tex
    check(f"pivot row {tex0}", f"{tex0} & {zp} & ${wp} \\pm {sw}$ & ${wt}$ & ${abs(dev) if round(dev, 2) == 0 else dev:.2f}\\sigma$ & ${float(g[15]):.2f}\\sigma$ & ${float(g[18]):.2f}\\sigma$")
piv = [abs(float(off[k][12])) for k in list(names)[:4]]
check("pivot agreement range", f"agrees to ${min(piv):.2f}$--${max(piv):.2f}\\sigma$")
bf = {k: off[k][6] for k in list(names)[:4]}
fmt_bf = []
for k in list(names)[:4]:
    lo, hi = map(float, off[k][6].split("-"))
    fmt_bf.append((lo, hi))
check("BF DESI+CMB", f"$B = {round(fmt_bf[0][0])}$--${round(fmt_bf[0][1])}$ (DESI+CMB)")
check("BF Pantheon+", f"${round(fmt_bf[1][0], -1):.0f}$--${round(fmt_bf[1][1], -1):.0f}$ (+Pantheon+)")
check("BF Union3", f"${round(fmt_bf[2][0])}$--${round(fmt_bf[2][1])}$ (+Union3)")
check("BF DESY5", f"${round(fmt_bf[3][0], -1):.0f}$--${round(fmt_bf[3][1], -1):.0f}$ (+DESY5)")
lbf = [off[k][7] for k in list(names)[:3]]
l0 = lbf[0].split("-"); l1 = lbf[1].split("-"); l2 = lbf[2].split("-")
check("BF LCDM:free", f"$B = {l0[0]}$--${l0[1]}$, ${float(l1[0]):.1f}$--${float(l1[1]):.1f}$ and ${l2[0]}$--${l2[1]}$")
ctx = R("v2_bf_context.txt")
pk = re.findall(r"^\S+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+) \|", ctx, re.M)
check("BF peak values", f"is ${round(float(pk[0][1]))}$, ${round(float(pk[1][1]))}$, ${round(float(pk[2][1]))}$ and ${round(float(pk[3][1]))}$")
ratios = [float(r[2]) for r in pk]
# B(phi:LCDM) = B(phi:free)/B(LCDM:free); conservative lower bound = min B(phi) / max B(LCDM), to 2 significant figures
lb = [fmt_bf[i][0] / float(off[k][7].split("-")[1]) for i, k in enumerate(list(names)[:3])]
lb = [f"{float(f'{v:.2g}'):.0f}" for v in lb]
check("BF phi:LCDM lower bounds", f"\\gtrsim {lb[0]}$, ${lb[1]}$ and ${lb[2]}$ for DESI+CMB, +Pantheon+ and +Union3")
ok = off[list(names)[3]][7].startswith("<") and re.search(r"^\+DESY5 \[baseline\].*not sampled", tb, re.M) is not None and "for DESY5 the \\LCDM\\ point lies outside the sampled region" in TEX
results.append(ok)
print(f"{'PASS' if ok else 'FAIL'}  DESY5 LCDM point outside the sampled region (stated, not 'zero')")
check("BF ratio to peak", f"reaches ${round(100 * min(ratios))}$--${round(100 * max(ratios))}\\%$")
fixed = re.search(r"Theory: w\(z=0.3\) = (-[\d.]+\d), w\(z=0.5\) = (-[\d.]+\d)", tb).groups()
check("fixed-z predictions (Table 3 caption)", f"$w(0.3) = {fixed[0]}$ and $w(0.5) = {fixed[1]}$")
chains = R("desi_chains_results.txt")
d3c = re.findall(r"Borse d=3 \(-1\+phi\^-3, -1/phi\)\s+[\d.]+%\s+([\d.]+)", chains)
check("Table 2 caption d=3 sigmas", f"lies at ${float(d3c[0]):.1f}$, ${float(d3c[1]):.1f}$, ${float(d3c[2]):.1f}$ and ${float(d3c[3]):.1f}\\sigma$")
nc = [float(x) for x in re.findall(r"non-crossing \(-phi/2, -0.191\)\s+[\d.]+%\s+[\d.]+\s+([\d.]+)", chains)]
check("non-crossing Gaussian range (official)", f"${rng(nc)}\\sigma$ (Gaussian estimate)")

# ---------------------------------------------------------------- background history (T12)
t12 = v2.split("T12.")[1].split("T11.")[0]
pk_z, lb = re.search(r"maximum at z = ([\d.]+) .*lookback time ([\d.]+) Gyr", t12).groups()
check("lookback time to z_cross", f"about ${lb}$\\,Gyr ago")
zeq, zacc = re.search(r"overtakes matter at z = ([\d.]+); acceleration begins at z = ([\d.]+)", t12).groups()
check("dark energy overtakes matter / acceleration onset", f"dark energy overtakes matter at $z = {zeq}$ and acceleration begins at $z = {zacc}$")
aend, aeq2, aw1 = re.search(r"ends at a = ([\d.]+); matter dominates again at a = ([\d.]+); w exceeds \+1 at a = ([\d.]+)", t12).groups()
check("future evolution numbers", f"universe is ${aend}$ times its present size, and matter would dominate again at ${float(aeq2):.1f}$ times")
check("CPL breakdown", f"beyond ${float(aw1):.1f}$ times the present size")

# ---------------------------------------------------------------- mu and growth extras
mu = R("mu_and_extras.txt")
m1 = re.search(r"\[gamma=1/phi\] mu-1 at z=0,0.5,1,2: \[(.*?)\]", mu).group(1).split(", ")
check("mu - 1 values (gamma = 1/phi)", f"$\\mu - 1 = {m1[0]}$, ${m1[1]}$, ${m1[2]}$ and ${m1[3]}$")
mu0 = float(re.search(r"\[gamma=1/phi\] mu0 = ([+-][\d.]+)", mu).group(1))
check("mu0 (gamma = 1/phi)", f"\\mu_0 \\simeq {mu0:.2f}")
res = float(re.search(r"\[gamma=1/phi\] mu0 = [+-][\d.]+; max residual = ([\d.]+)", mu).group(1))
check("mu fit residual", f"within ${res:.3f}$")
alt = [float(x.strip("'")) for x in re.search(r"predicted gamma \(GR\+F4\): \[(.*?)\]", mu).group(1).split(", ")]
check("alternative rule gamma values", f"$\\gamma = {alt[0]:.3f}, {alt[1]:.3f}, {alt[2]:.3f}$")
c = re.search(r"GR ([\d.]+)\s+1/phi ([\d.]+)\s+GR\+F4 ([\d.]+)", mu).groups()
check("chi2 1/phi vs GR+F4 (Remark)", f"$\\chi^2 = {c[1]}$ versus ${c[2]}$")
for f, row in (("growth_results.txt", "SDSS (6dF, MGS, BOSS, eBOSS)"), ("growth_results_desi.txt", "DESI DR1 (FS + PV)")):
    g = re.findall(r"Borse d=2\s+gamma free\s+[\d.]+\s+best gamma = ([\d.]+) -([\d.]+) \+([\d.]+)\s+\| 1/phi off by ([\d.]+)", R(f))
    (a, am, ap, asg), (b, bm, bp, bsg) = g
    check(f"growth table {row[:4]}", f"{row} & ${a}^{{+{ap}}}_{{-{am}}}$ & ${float(asg):.1f}\\sigma$ & ${b}^{{+{bp}}}_{{-{bm}}}$ & ${float(bsg):.1f}\\sigma$")
s8 = re.findall(r"Borse d=2\s+gamma = 1/phi\s+[\d.]+\s+[\d.]+\s+([\d.]+)", R("growth_results.txt"))[1]
check("implied S8", f"$S_8 = {float(s8):.3f}$")
gro = [float(x) for x in re.findall(r"Borse d=2\s+gamma free.*?1/phi off by ([\d.]+)", R("growth_results.txt") + R("growth_results_desi.txt"))][1::2]
check("growth CMB-anchored range (abstract)", f"growth index is ${rng(gro)}\\sigma$")

# ---------------------------------------------------------------- drag model
md = R("model_results.txt")
gc = dict(re.findall(r"xi =\s+([\d.]+): growth change \+([\d.]+)%", md))
check("drag: brief xi enhancement", f"$+{float(gc['0.332']):.2f}\\%$ for $\\xi = (7\\gr - 10)/4")
check("drag: xi=3 enhancement", f"$+{float(gc['3.000']):.1f}\\%$ for $\\xi = 3$")
tgt = abs(float(re.search(r"lowered by (-[\d.]+)%", md).group(1)))
check("drag: target suppression", f"The target is a {tgt:.1f}\\% lower growth factor")
ncc = [float(x) for x in re.findall(r"non-crossing boundary point disfavoured at ([\d.]+) sigma", md)]
check("non-crossing compressed range", f"excluded at ${rng(ncc)}\\sigma$ by the data of Table")
xi = float(re.search(r"xi needed for gamma = 1/phi suppression: ([\d.]+)", md).group(1))
bg = re.search(r"sigma_D/m_c ~ (\d+) b/GeV", md).group(1)
check("non-crossing xi and cross-section", f"$\\xi \\approx {xi:.1f}$ ($\\sigma_D/m_c \\approx {bg}$\\,b/GeV)")
dg = [float(x) for x in re.findall(r"z=[\d.]+: \+([\d.]+)", md.split("Q4.")[1].split("Delta gamma:")[1])]
check("non-crossing Delta gamma range", f"$\\Delta\\gamma = {min(dg):.3f}$--${max(dg):.3f}$")
tq4 = [float(x) for x in re.findall(r"\+([\d.]+)", md.split("Q4.")[1].split("target 1/phi - gamma_GR:")[1].split(")")[0])]
check("non-crossing required Delta gamma", f"close to the ${min(tq4):.3f}$--${max(tq4):.3f}$ required in that background")
tq2 = [float(x) for x in re.findall(r"z=[\d.]+: \+([\d.]+)", md.split("target Delta gamma = 1/phi - gamma_GR:")[1].split("\n")[0])]
check("crossing target gamma - gamma_GR range", f"$\\gamma - \\gamma_{{\\rm GR}} \\approx {min(tq2):.2f}$--${max(tq2):.2f}$ for $0 \\le z \\le 1.5$")

# ---------------------------------------------------------------- EDE, validation, CAMB, z*
r0 = R("results.txt")
rs, rd = re.search(r"r_s\(z\*\) = ([\d.]+) Mpc.*r_d = ([\d.]+) Mpc", r0).groups()
th = re.search(r"100 theta\* = ([\d.]+)", r0).group(1)
check("Planck validation", f"$r_s(z_*) = {rs}$\\,Mpc, $r_d = {rd}$\\,Mpc and $100\\theta_* = {th}$")
dr = [abs(float(x)) for x in re.findall(r"^\s+\d+\s+[\d.]+\s+(-[\d.]+)%", r0, re.M)]
check("EDE r_d reduction range", f"by only ${min(dr):.1f}$--${max(dr):.1f}\\%$")
fneed = float(re.search(r"r_d = 136.36 Mpc: f = ([\d.]+)", r0).group(1))
check("EDE fraction needed", f"$\\fede \\approx {fneed:.2f}$")
row = re.search(r"^\s+3500\s+[\d.]+\s+-[\d.]+%\s+([\d.]+)\s+([\d.]+)", r0, re.M).groups()
check("EDE H0 range", f"from ${float(row[1]):.1f}$ (holding $H_0 r_d$ fixed) to ${float(row[0]):.1f}$")
cb = R("camb_check.txt")
hd = max(abs(float(x)) for x in re.findall(r"^(?:H|D_M)\(z=[\d.]+\)\s+[\d.]+\s+[\d.]+\s+(-?[\d.]+)%", cb, re.M))
rr = abs(float(re.search(r"R \(shift parameter\)\s+[\d.]+\s+[\d.]+\s+(-?[\d.]+)%", cb).group(1)))
la = abs(float(re.search(r"l_A \(acoustic scale\)\s+[\d.]+\s+[\d.]+\s+(-?[\d.]+)%", cb).group(1)))
rdd = abs(float(re.search(r"r_d \(DESI fit vs CAMB\)\s+[\d.]+\s+[\d.]+\s+(-?[\d.]+)%", cb).group(1)))
gs = max(abs(float(x)) for x in re.findall(r"diff ([+-][\d.]+)%", cb))
check("CAMB H and D_M", f"agrees with CAMB to ${hd:.3f}\\%$ in $H(z)$ and $D_M(z)$")
check("CAMB R and l_A", f"${rr:.3f}\\%$ in $R$, ${la:.3f}\\%$ in $\\ell_A$")
check("CAMB r_d", f"${rdd:.2f}\\%$ in $r_d$")
check("CAMB growth shape", f"${gs:.2f}\\%$ in the shape of $f\\sigma_8(z)$")
zr = float(re.search(r"residual max ([\d.e-]+)", R("zstar_calibration.txt")).group(1))
ok = zr < 4e-7 and "residuals below $4\\times10^{-7}$" in TEX
results.append(ok)
print(f"{'PASS' if ok else 'FAIL'}  z* calibration residual < 4e-7 ({zr:.2e})")

# ---------------------------------------------------------------- joint EDE (E1a, E1b)
vb = R("v2b_results.txt")
e1a = vb.split("E1a.")[1].split("E1b.")[0]
rows = re.findall(r"^(DESI BAO|BAO\+\S+)\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)s\s+([\d.]+)s", e1a, re.M)
d1 = dict((r[0], r[1]) for r in rows)
check("E1a r_d-free phi vs free values", f"$\\gr$-ladder lies ${d1['DESI BAO']}$ (BAO only), ${d1['BAO+PantheonPlus']}$ (+Pantheon+), "
      f"${d1['BAO+Union3']}$ (+Union3), ${d1['BAO+DESY5']}$ (+DESY5) and ${d1['BAO+DESDovekie']}\\sigma$ (+DES-Dovekie)")
check("E1a LCDM range", f"\\LCDM\\ is disfavored at ${rng([float(r[2]) for r in rows])}\\sigma$. This agrees")
e1b = vb.split("E1b.")[1].split("B1.")[0]
rb = re.findall(r"^\s+\d+ DESI\+CMB\S*\s+([\d.]+)s\s+([\d.]+)s\s+([\d.]+)\s+([\d.]+)", e1b, re.M)
check("E1b phi vs free range", f"stays ${rng([float(r[0]) for r in rb])}\\sigma$ from the free fit")
check("E1b LCDM range", f"\\LCDM\\ is disfavored at ${rng([float(r[1]) for r in rb])}\\sigma$, and $H_0$ rises to ${rng([float(r[2]) for r in rb])}$")
shr = [100 * (1 - float(r[3])) for r in rb]
check("E1b r_d shrinkage range", f"($r_d$ shrinks by ${rng(shr)}\\%$)")

# ---------------------------------------------------------------- binned w(z) (B1 + calibration)
bc = R("binned_calibration.txt")
blab = {"DESI+CMB": "DESI+CMB", "+Pantheon+": "+Pantheon+", "+Union3": "+Union3", "+DESY5": "+DESY5"}
phis, cpls = [], []
for lab_ in blab:
    L = re.escape(lab_)
    ph = float(re.search(rf"^{L}\s+phi-ladder\s+[-\d., ]+?\s+[\d.]+\s+([\d.]+)\s*$", bc, re.M).group(1))
    cp = float(re.search(rf"^{L}\s+DESI CPL mean \([^)]*\)\s+[-\d., ]+?\s+[\d.]+\s+([\d.]+)\s*$", bc, re.M).group(1))
    phis.append(ph); cpls.append(cp)
    b1 = re.search(rf"^{re.escape(lab_ if lab_ != 'DESI+CMB' else 'DESI+CMB')}\s+bin0: (-[\d.]+)\+-([\d.]+).*?bin1: (-[\d.]+)\+-([\d.]+).*?bin2: (-[\d.]+)\+-([\d.]+)", vb, re.M).groups()
    check(f"binned table {lab_}", f"{lab_} & ${b1[0]} \\pm {b1[1]}$ & ${b1[2]} \\pm {b1[3]}$ & ${b1[4]} \\pm {b1[5]}$ & ${ph:.2f}\\sigma$ & ${cp:.2f}\\sigma$")
th = re.search(r"theory bin averages \(uniform in z\): \[0.0-0.7\] (-[\d.]+), \[0.7-1.4\] (-[\d.]+), \[1.4-2.1\] (-[\d.]+)", vb).groups()
check("binned theory averages", f"the predicted averages are ${th[0]}$, ${th[1]}$ and ${th[2]}$")
check("binned phi range", f"at ${min(phis):.1f}$--${max(phis):.2f}\\sigma$, comparable")
check("binned DESI CPL range", f"CPL models (${min(cpls):.1f}$--${max(cpls):.1f}\\sigma$)")

# ---------------------------------------------------------------- energy exchange (Open problem 4)
ee = R("energy_exchange_detail.txt")
xis = re.findall(r"xi = ([\d.]+) puts the crossing", ee)
check("energy-exchange couplings", f"with $\\xi = {float(xis[0]):.2f}$ and ${float(xis[1]):.2f}$ respectively")
z1 = re.findall(r"z=1.0: interacting w_eff = (-[\d.]+)\s+phi-ladder w = (-[\d.]+)", ee)
check("energy-exchange w_eff at z=1", f"(${float(z1[0][0]):.2f}$ and ${float(z1[1][0]):.2f}$ at $z = 1$, against ${float(z1[0][1]):.2f}$ predicted)")
negz = re.findall(r"non-positive at z = ([\d.]+)", ee)
check("energy-exchange negative density", f"turns negative by $z = {float(negz[0]):.1f}$ and ${float(negz[1]):.1f}$")
dev = []
for blk in ee.split("Q = ")[1:]:
    for z in ("0.0", "0.3", "0.5"):
        m = re.search(rf"z={z}: interacting w_eff = (-[\d.]+)\s+phi-ladder w = (-[\d.]+)", blk)
        dev.append(abs(float(m.group(1)) - float(m.group(2))))
ok = max(dev) <= 0.0405 and "to within $0.04$ up to $z = 0.5$" in TEX
results.append(ok)
print(f"{'PASS' if ok else 'FAIL'}  energy-exchange w_eff within 0.04 of prediction up to z=0.5 (max {max(dev):.3f})")

# ---------------------------------------------------------------- review additions (pass 4)
# DES-Dovekie published constraint (Popovic et al., arXiv:2511.07517 abstract; see VERIFICATION.md)
import math
PHI = (1 + 5 ** 0.5) / 2
dv = (abs(-PHI / 2 - (-0.803)) / 0.054, abs(-1 / PHI - (-0.72)) / 0.21)
check("prediction vs published DES-Dovekie (w0, wa)", f"the prediction lies ${dv[0]:.1f}\\sigma$ and ${dv[1]:.1f}\\sigma$ from these values")
check("robustness summary: binned", f"deviates by ${min(phis):.1f}$--${max(phis):.2f}\\sigma$, comparable to DESI's own best-fit CPL models (${min(cpls):.1f}$--${max(cpls):.1f}\\sigma$)")
check("robustness summary: no CMB, r_d free", f"the $\\gr$-ladder lies ${rng([float(v) for v in d1.values()])}\\sigma$ from the free fit; with CMB")
check("robustness summary: with EDE", f"it lies ${rng([float(r[0]) for r in rb])}\\sigma$ from it")
zps = [float(off[k][8]) for k in list(names)[:4]]
check("pivot redshift range", f"pivot redshifts (${min(zps):.2f}$--${max(zps):.2f}$)")
check("S8 in prediction section", f"$S_8 \\approx 0.81$ (${float(s8):.3f}$ in the best-fit background")
ab = TEX.split("\\begin{abstract}")[1].split("\\end{abstract}")[0].strip()
ok = len(ab) <= 1920
results.append(ok)
print(f"{'PASS' if ok else 'FAIL'}  abstract within arXiv's 1920-character limit ({len(ab)} characters of source)")

# ---------------------------------------------------------------- lock consistency
sha = open("PREDICTIONS.sha256").read().split()[0]
check("paper lock hash matches PREDICTIONS.sha256", "\\newcommand{\\lockhash}{" + sha + "}")

# ---------------------------------------------------------------- exact values (symbolic proofs)
sy = R("symbolic_proofs.txt")
w03, w05 = re.search(r"w\(0.3\) = .* = (-[\d.]+); w\(0.5\) = .* = (-[\d.]+)", sy).groups()
check("w(0.3) exact value", f"= {float(w03):.4f}")
check("w(0.5) exact value", f"= {float(w05):.4f}")
ok = re.search(r"(\d+)/(\d+) exact checks passed", sy)
ok = ok and ok.group(1) == ok.group(2)
results.append(bool(ok))
print(f"{'PASS' if ok else 'FAIL'}  all symbolic proofs pass")

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
