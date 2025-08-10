# Full Maxwell-Bloch Simulation for EIT-Based Quantum Optics

## Overview

This repository contains a comprehensive implementation of Maxwell-Bloch equations for simulating Electromagnetically Induced Transparency (EIT) in quantum optical systems. The project provides high-fidelity numerical simulations of light-matter interactions in atomic vapors, with particular focus on quantum memory applications.

## Features

- **Full Maxwell-Bloch Solver**: Implements coupled Maxwell-Bloch equations with Doppler averaging
- **Multiple Integration Methods**: Comparison between Euler and 4th-order Runge-Kutta methods
- **EIT Quantum Memory Simulations**: Storage and retrieval of quantum states
- **Co/Counter-propagating Geometries**: Analysis of different control field configurations
- **Retrieval Efficiency Analysis**: Heatmaps showing efficiency vs control parameters
- **High-Performance Computing**: Optimized for both CPU efficiency and accuracy

## Project Structure

```
├── src/                      # Source code
│   ├── core/                # Core physics modules
│   ├── solvers/            # Numerical solvers
│   └── visualization/      # Plotting utilities
├── examples/               # Example simulations
├── data/                   # Simulation results
├── figures/                # Generated figures
├── docs/                   # Documentation
└── tests/                  # Unit tests
```

## Installation

### Prerequisites
- Python 3.9 or higher
- NumPy, Pandas, Matplotlib, Pillow

### Setup

1. Clone the repository:
```bash
git clone https://github.com/alovladi007/Full-Maxwell-Bloch-Simulation-for-EIT-Based-Quantum-Optics.git
cd Full-Maxwell-Bloch-Simulation-for-EIT-Based-Quantum-Optics
```

2. Create a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Quick Start

Run all three main analyses:
```bash
python src/high_fidelity_qo_suite.py --outdir output --jpeg --euler_rk4 --heatmap --cocounter
```

## Usage Examples

### Basic Simulation
```python
from src.core import MBParams, maxwell_bloch_solver

# Set up parameters
params = MBParams(
    L=0.03,           # Medium length (m)
    Nz=120,           # Spatial points
    tmax=6e-6,        # Simulation time (s)
    Nt=4000,          # Time points
    N=5e16,           # Atomic density (m^-3)
    gamma_sg=2*pi*0.2e6  # Ground state dephasing
)

# Run simulation
results = maxwell_bloch_solver(params)
```

### High-Fidelity Analysis
```bash
python src/high_fidelity_qo_suite.py \
    --er_Nt 5000 --er_Nz 200 --er_Nv 15 \
    --hm_Nt 4000 --hm_Nz 180 --hm_doppler \
    --cc_Nt 4000 --cc_Nz 180 --cc_Nv 15
```

## Physics Background

The simulation solves the coupled Maxwell-Bloch equations for a three-level Λ-system:

- **Probe field**: Couples |g⟩ → |e⟩ transition
- **Control field**: Couples |s⟩ → |e⟩ transition
- **EIT window**: Created by destructive quantum interference

Key parameters:
- `γ_eg`: Excited state decay rate
- `γ_sg`: Ground state coherence decay rate
- `Ω_c`: Control field Rabi frequency
- `Δ_p`: Probe detuning

## Results

The simulation produces:
1. **Propagation dynamics**: Evolution of probe and control fields
2. **Storage/retrieval efficiency**: Quantum memory performance metrics
3. **Doppler effects**: Impact of atomic motion on EIT
4. **Comparison plots**: Different numerical methods and configurations

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this code in your research, please cite:
```bibtex
@software{maxwell_bloch_eit_2024,
  author = {Louis Antoine},
  title = {Full Maxwell-Bloch Simulation for EIT-Based Quantum Optics},
  year = {2024},
  url = {https://github.com/alovladi007/Full-Maxwell-Bloch-Simulation-for-EIT-Based-Quantum-Optics}
}
```

## Contact

For questions or collaborations, please contact: alovladi@gmail.com