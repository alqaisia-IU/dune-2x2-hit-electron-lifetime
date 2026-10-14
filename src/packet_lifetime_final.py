#!/usr/bin/env python3
"""
Preliminary 2x2 raw-packet / hit-based electron lifetime extraction.

For presentation use:
  raw packet file -> EXT-relative drift time -> pedestal-subtracted ADC
  -> near-threshold cut -> FWHM-restricted Landau-like x Gaussian MPV fits
  -> MPV(t) = A exp(-t/tau)

Recommended pedestal source:
  matching FLOW file, but ONLY to recover the static pedestal implicit in Q_raw:
      pedestal_ADC = raw_hit_ADC - Q_raw / ke_per_adc
No event building, reconstructed t0, Y-Z position, or corrected Q is used.

Important: nearline *.panel_ped.json is intentionally NOT used here because we
verified that it stores the same physics file's packet ADC mean/RMS, not a
separate baseline pedestal calibration.
"""

from pathlib import Path
from collections import defaultdict
import argparse, csv, json, math
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm

try:
    from scipy.stats import landau as LANDAU
    LANDAU_NAME = "scipy.stats.landau"
except Exception:
    from scipy.stats import moyal as LANDAU
    LANDAU_NAME = "scipy.stats.moyal (Landau approximation)"

KE_PER_ADC = 0.9414909055727555
TS_MOD = 2**31
MODULES = {"Module 0": (1,2), "Module 1": (3,4), "Module 2": (5,6), "Module 3": (7,8)}
GZ = np.linspace(-5, 5, 61)
GW = norm.pdf(GZ); GW /= GW.sum()


def ds(h5, path):
    obj = h5[path]
    if isinstance(obj, h5py.Dataset): return obj
    if isinstance(obj, h5py.Group) and "data" in obj: return obj["data"]
    raise RuntimeError(f"Cannot resolve dataset {path}")


def key(iog, ioc, chip, ch):
    return ((np.asarray(iog,dtype=np.uint64)<<24) |
            (np.asarray(ioc,dtype=np.uint64)<<16) |
            (np.asarray(chip,dtype=np.uint64)<<8) |
            np.asarray(ch,dtype=np.uint64))


def auto_flow(packet, version="v10"):
    s = str(packet)
    if "/nearline/packet/" not in s:
        raise RuntimeError("Cannot auto-derive FLOW path; pass --flow-file")
    s = s.replace("/nearline/packet/", f"/reflows/{version}/flow/", 1)
    if s.endswith(".h5"): s = s[:-3] + ".FLOW.hdf5"
    return Path(s)


def ext_counts(packet, ptype=7, chunk=2_000_000):
    out = defaultdict(int)
    with h5py.File(packet,"r") as f:
        p = f["packets"]
        for a in range(0,len(p),chunk):
            z = p.fields(["packet_type","io_group"])[a:min(a+chunk,len(p))]
            m = z["packet_type"] == ptype
            if np.any(m):
                u,n = np.unique(z["io_group"][m], return_counts=True)
                for x,y in zip(u,n): out[int(x)] += int(y)
    return dict(sorted(out.items()))


def extract(packet, ext_iog, tmax=200., tick=0.1, ext_ptype=7, chunk=2_000_000):
    """Match charge packets to the most recent EXT marker in FILE ORDER."""
    cols = ["packet_type","io_group","io_channel","chip_id","channel_id","timestamp","dataword"]
    out = defaultdict(list); last_ext = None
    with h5py.File(packet,"r") as f:
        p = f["packets"]; names = p.dtype.names or ()
        parity = "valid_parity" in names
        if parity: cols.append("valid_parity")
        for a in range(0,len(p),chunk):
            z = p.fields(cols)[a:min(a+chunk,len(p))]
            ei = np.flatnonzero((z["packet_type"]==ext_ptype) & (z["io_group"]==ext_iog))
            ets = z["timestamp"][ei].astype(np.int64)
            cm = z["packet_type"]==0
            if parity: cm &= z["valid_parity"].astype(bool)
            ci = np.flatnonzero(cm)
            if len(ci):
                j = np.searchsorted(ei, ci, side="right") - 1
                mt = np.full(len(ci), -1, dtype=np.int64)
                good = j>=0
                if np.any(good): mt[good] = ets[j[good]]
                if np.any(~good) and last_ext is not None: mt[~good] = last_ext
                ok = mt>=0
                ci2, mt2 = ci[ok], mt[ok]
                if len(ci2):
                    dt = np.mod(z["timestamp"][ci2].astype(np.int64)-mt2, TS_MOD)*tick
                    w = (dt>=0)&(dt<=tmax)
                    if np.any(w):
                        ii=ci2[w]; out["dt"].append(dt[w])
                        for c in ["io_group","io_channel","chip_id","channel_id","dataword"]:
                            out[c].append(z[c][ii].copy())
            if len(ets): last_ext=int(ets[-1])
    if not out["dt"]: raise RuntimeError("No packets in drift window; check EXT selection/timing")
    return {k:np.concatenate(v) for k,v in out.items()}


def flow_ped_map(flow, ke_per_adc=KE_PER_ADC, chunk=1_000_000):
    """Infer exact per-channel pedestal used by FLOW Q_raw."""
    sums=defaultdict(float); nums=defaultdict(int); mismatch=0
    with h5py.File(flow,"r") as f:
        r=ds(f,"charge/raw_hits"); c=ds(f,"charge/calib_prompt_hits")
        if len(r)!=len(c): raise RuntimeError("FLOW raw/calib hit lengths differ")
        for a in range(0,len(r),chunk):
            b=min(a+chunk,len(r))
            rr=r.fields(["id","ADC"])[a:b]
            cc=c.fields(["id","Q_raw","io_group","io_channel","chip_id","channel_id"])[a:b]
            same=rr["id"]==cc["id"]; mismatch += int(np.sum(~same))
            if not np.any(same): continue
            ped=rr["ADC"][same].astype(float)-cc["Q_raw"][same].astype(float)/ke_per_adc
            kk=key(cc["io_group"][same],cc["io_channel"][same],cc["chip_id"][same],cc["channel_id"][same])
            u,inv=np.unique(kk,return_inverse=True)
            ss=np.bincount(inv,weights=ped); nn=np.bincount(inv)
            for x,s,n in zip(u,ss,nn): sums[int(x)]+=float(s); nums[int(x)]+=int(n)
    if mismatch: raise RuntimeError(f"FLOW raw/calib IDs misaligned ({mismatch} rows)")
    return {k:sums[k]/nums[k] for k in sums}


def apply_ped(d,ped,min_q=6.,to_ke=False):
    kk=key(d["io_group"],d["io_channel"],d["chip_id"],d["channel_id"])
    pp=np.fromiter((ped.get(int(k),np.nan) for k in kk),float,count=len(kk))
    q=d["dataword"].astype(float)-pp
    m=np.isfinite(q)&(q>=min_q)
    scale=KE_PER_ADC if to_ke else 1.0
    return {"dt":d["dt"][m],"iog":d["io_group"][m],"q":q[m]*scale}, float(np.mean(np.isfinite(pp)))


def langau(x,A,loc,eta,sig):
    xx=np.asarray(x)[:,None]-sig*GZ[None,:]
    return A*(LANDAU.pdf(xx,loc=loc,scale=eta)@GW)


def fit_slice(q,nbins=80):
    q=np.asarray(q,float); q=q[np.isfinite(q)]
    if len(q)<300: return None
    lo,hi=np.percentile(q,[0.5,99.5])
    if hi<=lo: return None
    y,e=np.histogram(q,bins=nbins,range=(lo,hi)); x=(e[:-1]+e[1:])/2
    pk=int(np.argmax(y)); half=y[pk]/2
    L=pk; R=pk
    while L>0 and y[L]>=half: L-=1
    while R<len(y)-1 and y[R]>=half: R+=1
    L=max(0,min(L,pk-3)); R=min(len(y)-1,max(R,pk+3))
    xf=x[L:R+1]; yf=y[L:R+1].astype(float)
    if len(xf)<6 or yf.max()==0: return None
    width=max(xf[-1]-xf[0],1.)
    p0=[yf.max()*width,x[pk],max(width/5,.2),max(width/10,.1)]
    try:
        pop,cov=curve_fit(langau,xf,yf,p0=p0,
            bounds=([0,lo-width,.03,.03],[np.inf,hi+width,max(5*width,10),max(5*width,10)]),
            sigma=np.sqrt(np.maximum(yf,1)),absolute_sigma=True,maxfev=30000)
    except Exception: return None
    fine=np.linspace(lo,hi,1200); yy=langau(fine,*pop); mpv=float(fine[np.argmax(yy)])
    # Approximate MPV error from location covariance; enforce >= one fine-grid step.
    err=max(float(np.sqrt(max(cov[1,1],0))), float(fine[1]-fine[0]))
    chi=float(np.sum(((yf-langau(xf,*pop))/np.sqrt(np.maximum(yf,1)))**2)); ndf=max(len(xf)-4,1)
    return dict(mpv=mpv,err=err,n=len(q),chi2ndf=chi/ndf,pop=pop,hist=(y,e),fwhm=(L,R))


def expo(t,A,tau): return A*np.exp(-np.asarray(t)/tau)


def fit_tau(t,y,err,bootstrap=300,seed=123):
    t=np.asarray(t,float); y=np.asarray(y,float); err=np.asarray(err,float)
    m=np.isfinite(t)&np.isfinite(y)&np.isfinite(err)&(y>0)&(err>0); t,y,err=t[m],y[m],err[m]
    if len(t)<5: return None
    pop,cov=curve_fit(expo,t,y,p0=[y.max(),1500],bounds=([0,50],[np.inf,20000]),sigma=err,absolute_sigma=True,maxfev=30000)
    pred=expo(t,*pop); chi=float(np.sum(((y-pred)/err)**2)); ndf=max(len(t)-2,1); red=chi/ndf
    formal=float(np.sqrt(max(cov[1,1],0))); scaled=formal*math.sqrt(max(red,1))
    rng=np.random.default_rng(seed); resid=y-pred; vals=[]
    for _ in range(bootstrap):
        yb=pred+rng.choice(resid,size=len(resid),replace=True)
        try:
            pb,_=curve_fit(expo,t,yb,p0=pop,bounds=([0,50],[np.inf,20000]),maxfev=10000)
            vals.append(pb[1])
        except Exception: pass
    berr=float((np.percentile(vals,84)-np.percentile(vals,16))/2) if len(vals)>=20 else np.nan
    return dict(A=float(pop[0]),tau=float(pop[1]),formal_err=formal,scaled_err=scaled,bootstrap_err=berr,chi2ndf=red)


def analyze(dt,q,tmax=200.,bin_us=10.,nbins=80,bootstrap=300,seed=123):
    edges=np.arange(0,tmax+1e-9,bin_us)
    if len(edges)==0 or edges[-1] < tmax-1e-9: edges=np.append(edges,tmax)
    else: edges[-1]=tmax
    rows=[]; detail=[]
    for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        m=(dt>=a)&((dt<b) if i<len(edges)-2 else (dt<=b))
        r=fit_slice(q[m],nbins)
        rows.append(dict(t=(a+b)/2,lo=a,hi=b,n=int(m.sum()),mpv=np.nan if r is None else r["mpv"],err=np.nan if r is None else r["err"],slice_chi2ndf=np.nan if r is None else r["chi2ndf"]))
        detail.append(r)
    good=[r for r in rows if np.isfinite(r["mpv"]) and np.isfinite(r["err"])]
    life=fit_tau([r["t"] for r in good],[r["mpv"] for r in good],[r["err"] for r in good],bootstrap,seed) if len(good)>=5 else None
    return rows,detail,life


def life_text(l):
    if l is None: return "fit failed"
    e=l["bootstrap_err"] if np.isfinite(l["bootstrap_err"]) else l["scaled_err"]
    return rf"$\tau_e={l['tau']:.0f}\pm{e:.0f}\ \mu s$\n$\chi^2/ndf={l['chi2ndf']:.2f}$"


def plot_combined(out,title,d,rows,detail,life,unit):
    fig,ax=plt.subplots(1,2,figsize=(13,5))
    h=ax[0].hist2d(d["dt"],d["q"],bins=(80,100)); fig.colorbar(h[3],ax=ax[0],label="Packets")
    ax[0].set(xlabel=r"$t_{packet}-t_{EXT}$ [$\mu$s]",ylabel=f"Pedestal-subtracted charge [{unit}]",title="Selected raw charge packets")
    g=[r for r in rows if np.isfinite(r["mpv"])]
    if g:
        t=np.array([r["t"] for r in g]); y=np.array([r["mpv"] for r in g]); e=np.array([r["err"] for r in g])
        ax[1].errorbar(t,y,yerr=e,fmt="o",capsize=2,label="Langau MPV")
        if life:
            xx=np.linspace(0,t.max()*1.03,300); ax[1].plot(xx,expo(xx,life["A"],life["tau"]),label="Exponential fit")
            ax[1].text(.04,.96,life_text(life),transform=ax[1].transAxes,va="top")
    ax[1].set(xlabel=r"Drift time [$\mu$s]",ylabel=f"MPV [{unit}]",title="MPV vs drift time"); ax[1].grid(alpha=.25); ax[1].legend()
    fig.suptitle(title); fig.tight_layout(); fig.savefig(out/"lifetime_combined.png",dpi=180); plt.close(fig)


def plot_modules(out,title,res,unit):
    fig,axs=plt.subplots(2,2,figsize=(11,8)); axs=axs.ravel()
    for ax,(name,(rows,life)) in zip(axs,res.items()):
        g=[r for r in rows if np.isfinite(r["mpv"])]
        if g:
            t=np.array([r["t"] for r in g]); y=np.array([r["mpv"] for r in g]); e=np.array([r["err"] for r in g])
            ax.errorbar(t,y,yerr=e,fmt="o",capsize=2)
            if life:
                xx=np.linspace(0,t.max()*1.03,250); ax.plot(xx,expo(xx,life["A"],life["tau"])); ax.text(.04,.96,life_text(life),transform=ax.transAxes,va="top",fontsize=9)
        ax.set(title=name,xlabel=r"Drift time [$\mu$s]",ylabel=f"MPV [{unit}]"); ax.grid(alpha=.25)
    fig.suptitle(title+"\nPer-module raw-packet lifetime"); fig.tight_layout(); fig.savefig(out/"lifetime_modules.png",dpi=180); plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("packet_file",type=Path)
    ap.add_argument("--flow-file",default="auto")
    ap.add_argument("--reflow-version",default="v10")
    ap.add_argument("--outdir",type=Path,default=Path("packet_lifetime_final_out"))
    ap.add_argument("--ext-io-group",default="auto")
    ap.add_argument("--ext-packet-type",type=int,default=7)
    ap.add_argument("--tick-us",type=float,default=.1)
    ap.add_argument("--time-max-us",type=float,default=200.)
    ap.add_argument("--time-bin-us",type=float,default=10.)
    ap.add_argument("--min-q-adc",type=float,default=6.)
    ap.add_argument("--convert-to-ke",action="store_true")
    ap.add_argument("--hist-bins",type=int,default=80)
    ap.add_argument("--bootstrap",type=int,default=300)
    ap.add_argument("--chunk-size",type=int,default=2_000_000)
    a=ap.parse_args()

    packet=a.packet_file.resolve(); out=a.outdir.resolve(); out.mkdir(parents=True,exist_ok=True)
    counts=ext_counts(packet,a.ext_packet_type,a.chunk_size)
    if not counts: raise RuntimeError("No packet_type=7 EXT markers found")
    ext_iog=max(counts,key=counts.get) if a.ext_io_group=="auto" else int(a.ext_io_group)
    if counts.get(ext_iog,0)==0: raise RuntimeError(f"No EXT markers in io_group {ext_iog}; available {counts}")
    print("EXT counts by io_group:",counts); print("Using EXT io_group:",ext_iog)

    raw=extract(packet,ext_iog,a.time_max_us,a.tick_us,a.ext_packet_type,a.chunk_size)
    print("Packets in drift window:",len(raw["dt"])); print("Charge IO groups:",dict(zip(*np.unique(raw["io_group"],return_counts=True))))

    flow=auto_flow(packet,a.reflow_version) if a.flow_file=="auto" else Path(a.flow_file)
    if not flow.exists(): raise FileNotFoundError(f"FLOW pedestal source not found: {flow}")
    print("FLOW pedestal source:",flow)
    ped=flow_ped_map(flow,KE_PER_ADC,max(500_000,a.chunk_size//2)); print("Pedestal channels:",len(ped))
    d,cov=apply_ped(raw,ped,a.min_q_adc,a.convert_to_ke); print("Pedestal coverage:",cov); print("After threshold cut:",len(d["q"]))
    if len(d["q"])<5000: raise RuntimeError("Too few selected packets")

    unit="ke-" if a.convert_to_ke else "ADC"
    rows,detail,life=analyze(d["dt"],d["q"],a.time_max_us,a.time_bin_us,a.hist_bins,a.bootstrap,12345)
    if life: print(f"COMBINED tau = {life['tau']:.1f} us; bootstrap err = {life['bootstrap_err']:.1f} us; chi2/ndf = {life['chi2ndf']:.2f}")

    mods={}
    for name,groups in MODULES.items():
        m=np.isin(d["iog"],groups)
        if m.sum()<3000: continue
        rr,_,ll=analyze(d["dt"][m],d["q"][m],a.time_max_us,a.time_bin_us,a.hist_bins,a.bootstrap,12345+groups[0])
        mods[name]=(rr,ll)
        if ll: print(f"{name}: tau={ll['tau']:.1f} us  chi2/ndf={ll['chi2ndf']:.2f}")

    title=packet.stem+"\nPreliminary raw-packet/hit electron lifetime"
    plot_combined(out,title,d,rows,detail,life,unit); plot_modules(out,title,mods,unit)

    with open(out/"combined_mpv.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    summary={
      "packet_file":str(packet),"flow_pedestal_source":str(flow),"pedestal_coverage":cov,
      "ext_counts_by_io_group":counts,"ext_io_group":ext_iog,"tick_us":a.tick_us,
      "time_window_us":[0,a.time_max_us],"time_bin_us":a.time_bin_us,"min_q_adc":a.min_q_adc,
      "landau_implementation":LANDAU_NAME,
      "charge_convention":"raw packet ADC - FLOW-equivalent Q_raw pedestal",
      "method_status":"preliminary presentation-stage hit/packet method; no event building, no reconstructed t0, no Y-Z; FLOW used only to recover static Q_raw pedestal",
      "threshold_note":"min-q cut is an explicit approximation; exact historical Brooke/Russell single-trigger threshold-removal code has not been recovered",
      "combined_lifetime":life,
      "modules":{k:{"io_groups":list(MODULES[k]),"lifetime":v[1]} for k,v in mods.items()}
    }
    with open(out/"summary.json","w") as f: json.dump(summary,f,indent=2)
    print("Outputs:",out/"lifetime_combined.png",out/"lifetime_modules.png",out/"combined_mpv.csv",out/"summary.json",sep="\n  ")

if __name__=="__main__": main()
