# High-Fidelity Quantum Optics Analyses — Repro Kit

This kit lets you reproduce three analyses locally with adjustable fidelity:

1. **Euler vs RK4** Maxwell–Bloch propagation comparison (with Doppler)
2. **Retrieval efficiency heatmap** vs control ramp duration and ground-state dephasing \(\gamma_{sg}\)
3. **Co- vs counter-propagating** control across multiple Doppler widths

## Files
- `high_fidelity_qo_suite.py` — parameterized CLI script (Python 3.9+)
- `requirements.txt` — minimal dependencies
- This `README.md`

## Quick start
```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# Create an output folder named qo_outputs with all three analyses + JPEGs
python high_fidelity_qo_suite.py --outdir qo_outputs --jpeg --euler_rk4 --heatmap --cocounter
```

## Tunable parameters (CLI flags)
### Euler vs RK4
- `--er_Nt` (default 3000), `--er_Nz` (120), `--er_Nv` (9), `--er_sigma_v` (180.0 m/s)

### Retrieval heatmap
- `--hm_ramps_us "0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6"`
- `--hm_gammas_kHz "50,100,200,400,800"`
- `--hm_Nt` (default 2500), `--hm_Nz` (100)
- `--hm_doppler` to include Doppler (slower)

### Co vs Counter
- `--cc_sigmas "60,150,300"` (m/s)
- `--cc_Nt` (default 2500), `--cc_Nz` (100), `--cc_Nv` (9)

## Suggested presets
### Balanced quality (fast workstation)
```
--er_Nt 3000 --er_Nz 120 --er_Nv 9 --er_sigma_v 180 --hm_Nt 2500 --hm_Nz 100 --cc_Nt 2500 --cc_Nz 100 --cc_Nv 9
```

### High fidelity (can take many minutes)
```
--er_Nt 5000 --er_Nz 200 --er_Nv 15 --er_sigma_v 200 --hm_Nt 4000 --hm_Nz 180 --hm_doppler --cc_Nt 4000 --cc_Nz 180 --cc_Nv 15 --cc_sigmas "60,150,300,450"
```

## Outputs
Each analysis writes its own subfolder inside `--outdir` with:
- PNG figures (and JPEGs if `--jpeg` used)
- CSV data files for all plotted traces and matrices

## Notes & tips
- The integrator uses fixed-step RK4 with a mild limiter; for very stiff regimes (very large \(\gamma_{sg}\) or extreme detunings) consider smaller `dt` via larger `Nt`.
- Heatmap typically omits Doppler for speed; enable `--hm_doppler` if needed.
- All figures follow: **no seaborn, one chart per figure, default colors**.

Have fun exploring!