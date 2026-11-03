# Load-aware floating wake steering

Closed-loop wake steering of floating wind turbines (IEA 15 MW, VolturnUS-S) with FLORIDyn,
where the yaw setpoints are chosen by a cost function that trades off farm power against
tower and blade damage equivalent loads (DEL). FLORIDyn is coupled online to FAST.Farm
through ZeroMQ.

## Content

- `matlab/FLORIDyn/` - FLORIDyn code (MATLAB) extended to floating turbines: platform pitch and
  vertical wake deflection, closed-loop controller with the load-aware cost function, EnKF,
  DEL surrogate models (FFNN, optional GPR correction).
  - `Simulations/2026_3T_DEL/3T_IEA15_floating_00/`: three-turbine case study (settings,
    `CLC_init_floating_fun_ZMQ_<run>.m`, `postprocess/` Python scripts for power and DEL figures).
- `controller_floris/` - Python side of the ZeroMQ link with FAST.Farm (`FLORIDyn_<run>.py`).
- `fastfarm_simulation/` - FAST.Farm inputs of two example cases: run 15001 (FFNN+FFNN) and run 15013
  (same setup, tower second stage = GPR instead of FFNN, `residualTowerGPR_wd.m`): `<run>.fstf`,
  turbine inputs in `wind_turbines/T1..T3`, inflow file, MLife settings and `.fstf` generators.
  Only the folder structure and the inputs are included, not the outputs.
- `scripts/` - PowerShell launchers of the run families.

## Not included

Simulation outputs (`.out`, `.outb`, `Results/`, MLife workbooks), the trained GPR residual models
(`gpr_residual_*.mat`, larger than GitHub's limit) and the FAST.Farm / OpenFAST executables.
Paths inside the scripts are absolute paths of the author's machines and must be adapted.

The compiled ROSCO controller libraries used by the three turbines (`libdiscon.T*.so`, Linux, in
`wind_turbines/T*/turbine_controllers/ins/`) are included; ROSCO is distributed under its own license
(Apache 2.0, <https://github.com/NREL/ROSCO>).

## Requirements

MATLAB, Python 3 (numpy, pandas, matplotlib, openpyxl) for `postprocess/`, FAST.Farm with ZeroMQ
support and ROSCO.
