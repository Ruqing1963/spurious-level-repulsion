"""
spurious_goe.py
---------------
Three ways finite/structured data fabricate GOE-like level repulsion in the
consecutive spacing ratio <r>, with the correct controls. Reproduces all numbers
and figures for the note. No third-party data download required: the SSW central
dates (public, NOAA SSW Compendium) are embedded; everything else is synthetic.
"""
import numpy as np, datetime as dt, collections, csv, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import poisson

R_POISSON, R_GOE, R_GUE = 0.386, 0.531, 0.603

def r_of_spacings(s):
    s = np.asarray(s, float); s = s[s > 0]
    return float((np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])).mean())

def r_of_points(x):
    return r_of_spacings(np.diff(np.sort(x)))

# ---------- Hazard 1: minimum-separation declustering ----------
def declustering_curve(n=400, fracs=(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7), seed=0):
    rng = np.random.default_rng(seed); out = []
    for f in fracs:
        rr = []
        for _ in range(200):
            x = np.sort(rng.uniform(0, n, n))           # Poisson points
            med = np.median(np.diff(x)); keep = [x[0]]
            for xi in x[1:]:
                if xi - keep[-1] >= f*med: keep.append(xi)  # greedy min-sep declustering
            if len(keep) > 5: rr.append(r_of_points(np.array(keep)))
        out.append((f, np.mean(rr)))
    return out

# ---------- Hazard 2: seasonal gating (SSW) ----------
SSW = ["1958-01-30","1960-01-17","1963-01-30","1965-12-18","1966-02-23","1968-01-07",
"1968-11-29","1970-01-02","1971-01-18","1971-03-20","1973-01-31","1977-01-09","1979-02-22",
"1980-02-29","1981-02-06","1981-03-04","1981-12-04","1984-02-24","1985-01-01","1987-01-23",
"1987-12-08","1988-03-14","1989-02-21","1998-12-15","1999-02-26","2000-03-20","2001-02-11",
"2001-12-31","2003-01-18","2004-01-05","2006-01-21","2007-02-24","2008-02-22","2009-01-24",
"2010-02-09","2010-03-24","2013-01-07","2018-02-12","2019-01-02","2021-01-05","2023-02-16"]
def ssw_analysis(seed=0):
    rng = np.random.default_rng(seed)
    d = [dt.date.fromisoformat(s) for s in SSW]
    o = np.array([x.toordinal() for x in d], float); s = np.diff(o)
    ro = r_of_spacings(s)
    rs = np.mean([r_of_spacings(rng.permutation(s)) for _ in range(3000)])
    win = [(x.year if x.month >= 7 else x.year-1) for x in d]
    cnt = collections.Counter(win); yrs = range(min(win), max(win)+1)
    counts = np.array([cnt.get(y, 0) for y in yrs])
    vm = counts.var()/counts.mean()
    boot = [np.var(c)/np.mean(c) for c in (rng.choice(counts, len(counts), replace=True) for _ in range(5000))]
    return dict(intervals=s, r_obs=ro, r_shuf=rs, counts=counts, mean=counts.mean(),
                vm=vm, vm_ci=(np.percentile(boot,2.5), np.percentile(boot,97.5)),
                frac_short=(s<90).mean(), frac_long=(s>270).mean())

# ---------- Hazard 3: Wishart / sample-correlation null ----------
def wishart_r(P, T, trim=0.15, reps=40, seed=0):
    rng = np.random.default_rng(seed); rr = []
    for _ in range(reps):
        X = rng.standard_normal((P, T))             # P UNCORRELATED series, length T
        ev = np.linalg.eigvalsh(np.corrcoef(X))
        k = int(len(ev)*trim); rr.append(r_of_points(ev[k:len(ev)-k]))
    return float(np.mean(rr))

def main(outdir):
    os.makedirs(f"{outdir}/results", exist_ok=True); os.makedirs(f"{outdir}/figures", exist_ok=True)
    dec = declustering_curve(); ssw = ssw_analysis()
    wish = [("P=100,T=400", wishart_r(100,400)), ("P=100,T=200", wishart_r(100,200)),
            ("P=150,T=600", wishart_r(150,600)), ("P=200,T=300", wishart_r(200,300))]
    # results CSV
    with open(f"{outdir}/results/summary.csv","w",newline="") as f:
        w=csv.writer(f); w.writerow(["hazard","quantity","value"])
        for fr,rv in dec: w.writerow(["declustering",f"r at min_sep_frac={fr}",round(rv,3)])
        w.writerow(["seasonal_gating","SSW raw interval <r>_obs",round(ssw['r_obs'],3)])
        w.writerow(["seasonal_gating","SSW raw interval <r>_shuf",round(ssw['r_shuf'],3)])
        w.writerow(["seasonal_gating","frac intervals <90d",round(ssw['frac_short'],3)])
        w.writerow(["seasonal_gating","frac intervals >270d",round(ssw['frac_long'],3)])
        w.writerow(["seasonal_gating","de-seasonalized var/mean",round(ssw['vm'],3)])
        w.writerow(["seasonal_gating","var/mean CI low",round(ssw['vm_ci'][0],3)])
        w.writerow(["seasonal_gating","var/mean CI high",round(ssw['vm_ci'][1],3)])
        for lab,rv in wish: w.writerow(["wishart_null",f"iid bulk <r> {lab}",round(rv,3)])
    print("declustering:", [(f,round(r,3)) for f,r in dec])
    print(f"SSW: raw <r>_obs={ssw['r_obs']:.3f} shuf={ssw['r_shuf']:.3f}; var/mean={ssw['vm']:.3f} CI={tuple(round(x,2) for x in ssw['vm_ci'])}")
    print("Wishart iid bulk <r>:", [(l,round(r,3)) for l,r in wish])

    # ---- Figure 1: three hazards ----
    fig,ax=plt.subplots(1,3,figsize=(15,4.4))
    fr=[f for f,_ in dec]; rv=[r for _,r in dec]
    ax[0].plot(fr,rv,'o-',color='#C0392B'); ax[0].axhline(R_GOE,ls='--',color='#2471A3',label='GOE 0.531')
    ax[0].axhline(R_POISSON,ls='--',color='#888',label='Poisson 0.386')
    ax[0].set_xlabel('declustering strength (min-sep / median)'); ax[0].set_ylabel(r'$\langle r\rangle$')
    ax[0].set_title('(a) Hazard 1: declustering fabricates GOE\n(input = pure random points)'); ax[0].legend(fontsize=8)
    ax[1].hist(ssw['intervals'],bins=np.logspace(1.3,3.6,20),color='#7D3C98',alpha=.8)
    ax[1].set_xscale('log'); ax[1].set_xlabel('SSW inter-event interval (days)'); ax[1].set_ylabel('count')
    ax[1].axvline(365,ls=':',color='k'); ax[1].set_title(f'(b) Hazard 2: seasonal gating\nbimodal intervals; naive <r>={ssw["r_obs"]:.2f} (=shuf, artifact)')
    labs=[l for l,_ in wish]; rws=[r for _,r in wish]
    ax[2].bar(range(len(labs)),rws,color='#27AE60'); ax[2].set_xticks(range(len(labs))); ax[2].set_xticklabels(labs,rotation=30,fontsize=7,ha='right')
    ax[2].axhline(R_GOE,ls='--',color='#2471A3',label='GOE 0.531'); ax[2].axhline(R_POISSON,ls='--',color='#888',label='Poisson 0.386')
    ax[2].set_ylim(0.35,0.58); ax[2].set_ylabel(r'bulk $\langle r\rangle$')
    ax[2].set_title('(c) Hazard 3: Wishart null\n(iid data -> sample corr. matrix -> GOE)'); ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_hazards.pdf",bbox_inches='tight'); fig.savefig(f"{outdir}/figures/fig_hazards.png",dpi=110,bbox_inches='tight')

    # ---- Figure 2: SSW de-seasonalized count vs Poisson ----
    fig,axn=plt.subplots(1,1,figsize=(6.4,4.6))
    c=ssw['counts']; nwin=len(c); m=ssw['mean']; ks=np.arange(0,4)
    obs=[ (c==k).sum() for k in ks ]; exp=[ nwin*poisson.pmf(k,m) for k in ks ]
    w=0.38
    axn.bar(ks-w/2,obs,w,label='observed',color='#2980B9'); axn.bar(ks+w/2,exp,w,label=f'Poisson({m:.2f})',color='#bbb')
    axn.set_xlabel('major SSWs per winter'); axn.set_ylabel('number of winters'); axn.set_xticks(ks)
    axn.set_title(f'De-seasonalized SSW counts: under-dispersed (var/mean={ssw["vm"]:.2f},\n95% CI {tuple(round(float(x),2) for x in ssw["vm_ci"])} excludes 1) = weak refractory clock, not GOE')
    axn.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_ssw_counts.pdf",bbox_inches='tight'); fig.savefig(f"{outdir}/figures/fig_ssw_counts.png",dpi=110,bbox_inches='tight')
    print("figures written")

if __name__=="__main__":
    import argparse; ap=argparse.ArgumentParser(); ap.add_argument("--outdir",default=".."); a=ap.parse_args()
    main(a.outdir)
