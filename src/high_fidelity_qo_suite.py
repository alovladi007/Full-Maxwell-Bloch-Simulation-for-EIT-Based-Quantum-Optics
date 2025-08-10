#!/usr/bin/env python3
"""
High-Fidelity Quantum Optics Analyses Suite
===========================================

This script reproduces three analyses with configurable fidelity:

1) Euler vs RK4 Maxwell–Bloch propagation comparison (with Doppler)
2) Retrieval efficiency heatmap vs control ramp duration and gamma_sg
3) Co- vs counter-propagating control at multiple Doppler widths

Notes
-----
- Uses matplotlib only (no seaborn), one figure per plot, no style overrides.
- Figures are saved as PNG; optional JPEG compression is available.
- Designed for local execution; high-fidelity settings may take minutes.
"""

import os
import argparse
from math import pi, sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

# ------------------------ Core Physics Utilities ------------------------

c0 = 299792458.0
eps0 = 8.854e-12
hbar = 1.054e-34

class MBParams:
    def __init__(self, L=0.03, Nz=120, tmax=6e-6, Nt=4000,
                 lam_p=780e-9, lam_c=780e-9,
                 N=5e16, mu=3.0e-29,
                 gamma_eg=2*pi*6e6, gamma_sg=2*pi*0.2e6,
                 Delta_p=0.0, delta_c=0.0,
                 sigma_v=180.0, Nv=11, control_geom="co"):
        self.L=L; self.Nz=Nz; self.tmax=tmax; self.Nt=Nt
        self.lam_p=lam_p; self.lam_c=lam_c
        self.N=N; self.mu=mu
        self.gamma_eg=gamma_eg; self.gamma_sg=gamma_sg
        self.Delta_p=Delta_p; self.delta_c=delta_c
        self.sigma_v=sigma_v; self.Nv=Nv; self.control_geom=control_geom
        self.kp = 2*pi/lam_p; self.kc = 2*pi/lam_c
        self.eta = (self.N * self.mu**2 * (2*pi*c0/self.lam_p)) / (2*eps0*hbar*c0)

def gaussian_pulse(t, t0, tau_fwhm, area=0.5 * 2*pi*1e6):
    sigma = tau_fwhm / (2.0*sqrt(2*np.log(2.0)))
    env = np.exp(-0.5*((t - t0)/sigma)**2)
    A = np.trapz(env, t)
    return (area / A) * env

def control_const(t, Omega_c0):
    return np.full_like(t, Omega_c0)

def control_ramp_storage(t, Omega_c0=2*pi*10e6, t0_on=2.0e-6, ramp=0.8e-6, hold=0.8e-6):
    Oc = np.zeros_like(t)
    t1=t0_on; t2=t1+ramp; t3=t2+hold; t4=t3+ramp
    for i,ti in enumerate(t):
        if ti < t1: Oc[i]=Omega_c0
        elif ti < t2: Oc[i]=Omega_c0*(1-(ti-t1)/ramp)
        elif ti < t3: Oc[i]=0.0
        elif ti < t4: Oc[i]=Omega_c0*((ti-t3)/ramp)
        else: Oc[i]=Omega_c0
    return Oc

def doppler_weights(Nv, sigma_v):
    vs = np.linspace(-3*sigma_v, 3*sigma_v, Nv)
    ws = np.exp(-0.5*(vs/sigma_v)**2)
    ws /= ws.sum()
    return vs, ws

def atomic_rhs(y, Op, Oc, Dp, dtwo, g_eg, g_sg):
    rho_eg = y[0] + 1j*y[1]
    rho_sg = y[2] + 1j*y[3]
    drho_eg = -(g_eg/2 + 1j*Dp)*rho_eg + 1j*Op/2 + 1j*Oc/2 * rho_sg
    drho_sg = -(g_sg + 1j*dtwo)*rho_sg + 1j*Oc/2 * rho_eg
    return np.array([drho_eg.real, drho_eg.imag, drho_sg.real, drho_sg.imag])

def rk4_integrate(t, Op_t, Oc_t, Dp, dtwo, g_eg, g_sg, limiter=5.0):
    y = np.zeros(4)
    dt = t[1]-t[0]
    rho_eg = np.zeros_like(t, dtype=np.complex128)
    for i in range(len(t)-1):
        Op=Op_t[i]; Oc=Oc_t[i]
        k1=atomic_rhs(y,Op,Oc,Dp,dtwo,g_eg,g_sg)
        k2=atomic_rhs(y+0.5*dt*k1,Op,Oc,Dp,dtwo,g_eg,g_sg)
        k3=atomic_rhs(y+0.5*dt*k2,Op,Oc,Dp,dtwo,g_eg,g_sg)
        k4=atomic_rhs(y+dt*k3,Op,Oc,Dp,dtwo,g_eg,g_sg)
        y=y+(dt/6.0)*(k1+2*k2+2*k3+k4)
        y=np.clip(y,-limiter,limiter)
        rho_eg[i+1]=y[0]+1j*y[1]
    return rho_eg

def mb_propagate_rk4(params, Op_in, Oc, include_doppler=True):
    Nt=params.Nt; Nz=params.Nz
    t = np.linspace(0, params.tmax, Nt)
    z = np.linspace(0, params.L, Nz+1)
    dz = z[1]-z[0]
    vs, ws = doppler_weights(params.Nv, params.sigma_v) if include_doppler else (np.array([0.0]), np.array([1.0]))
    k2ph = params.kp - params.kc if params.control_geom=="co" else params.kp + params.kc
    Op = np.zeros((Nz+1, Nt), dtype=np.complex128)
    Op[0,:] = Op_in.astype(np.complex128)
    for iz in range(Nz):
        rho_eg_avg = np.zeros(Nt, dtype=np.complex128)
        for v, w in zip(vs, ws):
            Dp = params.Delta_p - params.kp * v
            dtwo = (params.Delta_p - params.delta_c) - k2ph * v
            rho_eg = rk4_integrate(t, Op[iz,:], Oc, Dp, dtwo, params.gamma_eg, params.gamma_sg)
            rho_eg_avg += w * rho_eg
        Op[iz+1,:] = Op[iz,:] + 1j * params.eta * rho_eg_avg * dz
    return t, z, Op

def mb_propagate_euler(params, Op_in, Oc, include_doppler=True):
    Nt=params.Nt; Nz=params.Nz
    t = np.linspace(0, params.tmax, Nt)
    dt = t[1]-t[0]
    z = np.linspace(0, params.L, Nz+1)
    dz = z[1]-z[0]
    vs, ws = doppler_weights(params.Nv, params.sigma_v) if include_doppler else (np.array([0.0]), np.array([1.0]))
    k2ph = params.kp - params.kc if params.control_geom=="co" else params.kp + params.kc
    Op = np.zeros((Nz+1, Nt), dtype=np.complex128)
    Op[0,:] = Op_in.astype(np.complex128)
    for iz in range(Nz):
        rho_eg_avg = np.zeros(Nt, dtype=np.complex128)
        rho_sg = np.zeros(Nt, dtype=np.complex128)
        for v, w in zip(vs, ws):
            Dp = params.Delta_p - params.kp * v
            dtwo = (params.Delta_p - params.delta_c) - k2ph * v
            rho_eg = np.zeros(Nt, dtype=np.complex128)
            rho_sg = np.zeros(Nt, dtype=np.complex128)
            for i in range(Nt-1):
                drho_eg = (-(params.gamma_eg/2 + 1j*Dp)*rho_eg[i] + 1j*Op[iz,i]/2 + 1j*Oc[i]/2 * rho_sg[i]) * dt
                drho_sg = (-(params.gamma_sg + 1j*dtwo)*rho_sg[i] + 1j*Oc[i]/2 * rho_eg[i]) * dt
                rho_eg[i+1] = rho_eg[i] + drho_eg
                rho_sg[i+1] = rho_sg[i] + drho_sg
            rho_eg_avg += w * rho_eg
        Op[iz+1,:] = Op[iz,:] + 1j * params.eta * rho_eg_avg * dz
    return t, z, Op

# ------------------------ Figure & IO Utilities ------------------------

def save_png(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def compress_to_jpeg(png_path, jpeg_path, quality=80, max_width=1800):
    img = Image.open(png_path)
    w, h = img.size
    if w > max_width:
        scale = max_width / w
        img = img.resize((int(w*scale), int(h*scale)), Image.Resampling.LANCZOS)
    img.convert("RGB").save(jpeg_path, "JPEG", quality=quality)

# ------------------------ Analyses ------------------------

def analysis_euler_vs_rk4(outdir, Nt=3000, Nz=120, Nv=9, sigma_v=180.0):
    os.makedirs(outdir, exist_ok=True)
    p = MBParams(Nt=Nt, Nz=Nz, Nv=Nv, sigma_v=sigma_v)
    t = np.linspace(0, p.tmax, p.Nt)
    Op_in = gaussian_pulse(t, 2.2e-6, 0.8e-6)
    Oc = control_const(t, 2*pi*8e6)
    t_e, _, Op_e = mb_propagate_euler(p, Op_in, Oc, include_doppler=True)
    t_r, _, Op_r = mb_propagate_rk4(p, Op_in, Oc, include_doppler=True)

    # Figure
    fig = plt.figure()
    plt.plot(t*1e6, np.abs(Op_e[-1,:])**2, label="Euler")
    plt.plot(t*1e6, np.abs(Op_r[-1,:])**2, label="RK4")
    plt.xlabel("Time (µs)"); plt.ylabel(r"Probe intensity $\propto |\Omega_p|^2$")
    plt.title("MB Propagation: Euler vs RK4 (high fidelity)"); plt.legend()
    png = os.path.join(outdir, "Fig_Euler_vs_RK4.png")
    save_png(fig, png)

    # Data
    df = pd.DataFrame({"time_us": t*1e6,
                       "output_euler_abs2": np.abs(Op_e[-1,:])**2,
                       "output_rk4_abs2": np.abs(Op_r[-1,:])**2})
    df.to_csv(os.path.join(outdir, "mb_compare_euler_rk4.csv"), index=False)

def analysis_retrieval_heatmap(outdir, ramps_us, gammas_kHz, Nt=2500, Nz=100, include_doppler=False):
    os.makedirs(outdir, exist_ok=True)
    p = MBParams(Nt=Nt, Nz=Nz)
    t = np.linspace(0, p.tmax, p.Nt)
    Op_in = gaussian_pulse(t, 2.2e-6, 0.8e-6)
    Ein = np.trapz(np.abs(Op_in)**2, t)

    eff = np.zeros((len(gammas_kHz), len(ramps_us)))
    for i, gk in enumerate(gammas_kHz):
        p.gamma_sg = 2*pi*gk*1e3
        for j, ru in enumerate(ramps_us):
            Oc = control_ramp_storage(t, Omega_c0=2*pi*10e6, t0_on=2.0e-6, ramp=ru*1e-6, hold=0.9e-6)
            _,_,Op = mb_propagate_rk4(p, Op_in, Oc, include_doppler=include_doppler)
            t_retrieve_start = 2.0e-6 + ru*1e-6 + 0.9e-6
            idx = np.argmin(np.abs(t - t_retrieve_start))
            Eout = np.trapz(np.abs(Op[-1, idx:])**2, t[idx:])
            eff[i, j] = Eout / Ein

    # Save CSV
    df = pd.DataFrame(eff, index=[f"{g:.0f}kHz" for g in gammas_kHz],
                      columns=[f"{ru:.2f}" for ru in ramps_us])
    df.to_csv(os.path.join(outdir, "retrieval_efficiency_heatmap.csv"))

    # Figure
    fig = plt.figure()
    plt.imshow(eff, origin="lower", aspect="auto",
               extent=[ramps_us[0], ramps_us[-1], gammas_kHz[0], gammas_kHz[-1]])
    plt.colorbar(label="Retrieval efficiency")
    plt.xlabel("Ramp duration (µs)"); plt.ylabel(r"$\gamma_{sg}/2\pi$ (kHz)")
    plt.title("EIT Retrieval Efficiency vs Ramp & Dephasing")
    png = os.path.join(outdir, "Fig_Retrieval_Heatmap.png")
    save_png(fig, png)

def analysis_co_vs_counter_multi_sigma(outdir, sigmas, Nt=2500, Nz=100, Nv=9):
    os.makedirs(outdir, exist_ok=True)
    curves = []
    for s in sigmas:
        p_co = MBParams(Nt=Nt, Nz=Nz, Nv=Nv, sigma_v=s, control_geom="co")
        p_ct = MBParams(Nt=Nt, Nz=Nz, Nv=Nv, sigma_v=s, control_geom="counter")
        t = np.linspace(0, p_co.tmax, p_co.Nt)
        Op_in = gaussian_pulse(t, 2.2e-6, 0.8e-6)
        Oc = control_const(t, 2*pi*8e6)
        _,_,Op_co = mb_propagate_rk4(p_co, Op_in, Oc, include_doppler=True)
        _,_,Op_ct = mb_propagate_rk4(p_ct, Op_in, Oc, include_doppler=True)
        curves.append((s, t, Op_co[-1,:], Op_ct[-1,:]))
        pd.DataFrame({"time_us": t*1e6,
                      "output_co_abs2": np.abs(Op_co[-1,:])**2,
                      "output_counter_abs2": np.abs(Op_ct[-1,:])**2}).to_csv(
            os.path.join(outdir, f"co_counter_sigma_{int(s)}ms.csv"), index=False)

    # Figure
    fig = plt.figure()
    for s, t, co, ct in curves:
        plt.plot(t*1e6, np.abs(co)**2, label=f"Co, σv={int(s)} m/s")
        plt.plot(t*1e6, np.abs(ct)**2, linestyle="--", label=f"Counter, σv={int(s)} m/s")
    plt.xlabel("Time (µs)"); plt.ylabel(r"Probe intensity $\propto |\Omega_p|^2$")
    plt.title("Co vs Counter at Multiple Doppler Widths")
    plt.legend(ncol=2, fontsize=9)
    png = os.path.join(outdir, "Fig_Co_vs_Counter_multi_sigma.png")
    save_png(fig, png)

# ------------------------ CLI ------------------------

def main():
    ap = argparse.ArgumentParser(description="High-Fidelity Quantum Optics Analyses Suite")
    ap.add_argument("--outdir", type=str, default="qo_outputs", help="Output directory")
    ap.add_argument("--jpeg", action="store_true", help="Also save JPEG versions of figures")
    # Euler vs RK4
    ap.add_argument("--euler_rk4", action="store_true", help="Run Euler vs RK4 comparison")
    ap.add_argument("--er_Nt", type=int, default=3000)
    ap.add_argument("--er_Nz", type=int, default=120)
    ap.add_argument("--er_Nv", type=int, default=9)
    ap.add_argument("--er_sigma_v", type=float, default=180.0)
    # Retrieval heatmap
    ap.add_argument("--heatmap", action="store_true", help="Run retrieval efficiency heatmap")
    ap.add_argument("--hm_ramps_us", type=str, default="0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6")
    ap.add_argument("--hm_gammas_kHz", type=str, default="50,100,200,400,800")
    ap.add_argument("--hm_Nt", type=int, default=2500)
    ap.add_argument("--hm_Nz", type=int, default=100)
    ap.add_argument("--hm_doppler", action="store_true", help="Include Doppler in heatmap (slower)")
    # Co vs Counter
    ap.add_argument("--cocounter", action="store_true", help="Run co vs counter at multiple Doppler widths")
    ap.add_argument("--cc_sigmas", type=str, default="60,150,300")
    ap.add_argument("--cc_Nt", type=int, default=2500)
    ap.add_argument("--cc_Nz", type=int, default=100)
    ap.add_argument("--cc_Nv", type=int, default=9)

    args = ap.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Run analyses based on flags
    if args.euler_rk4:
        analysis_euler_vs_rk4(outdir=os.path.join(outdir, "euler_vs_rk4"),
                              Nt=args.er_Nt, Nz=args.er_Nz, Nv=args.er_Nv, sigma_v=args.er_sigma_v)

    if args.heatmap:
        ramps_us = [float(x) for x in args.hm_ramps_us.split(",") if x.strip()]
        gammas_kHz = [float(x) for x in args.hm_gammas_kHz.split(",") if x.strip()]
        analysis_retrieval_heatmap(outdir=os.path.join(outdir, "retrieval_heatmap"),
                                   ramps_us=ramps_us, gammas_kHz=gammas_kHz,
                                   Nt=args.hm_Nt, Nz=args.hm_Nz, include_doppler=args.hm_doppler)

    if args.cocounter:
        sigmas = [float(x) for x in args.cc_sigmas.split(",") if x.strip()]
        analysis_co_vs_counter_multi_sigma(outdir=os.path.join(outdir, "co_vs_counter"),
                                           sigmas=sigmas, Nt=args.cc_Nt, Nz=args.cc_Nz, Nv=args.cc_Nv)

    # Optional JPEGs
    if args.jpeg:
        for root, _, files in os.walk(outdir):
            for fn in files:
                if fn.lower().endswith(".png"):
                    png_path = os.path.join(root, fn)
                    jpg_path = os.path.splitext(png_path)[0] + ".jpg"
                    compress_to_jpeg(png_path, jpg_path)

    print("Done. Outputs in:", outdir)

if __name__ == "__main__":
    main()
