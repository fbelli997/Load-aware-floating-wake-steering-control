import os
import re
import numpy as np
from rosco.toolbox.control_interface import wfc_zmq_server
import pandas as pd
from pathlib import Path
from scipy.signal import butter, filtfilt
import json
import time as time_mod

print('********************CIAO: SONO IL CONTROLLER CHE CHIAMA FLORIDyn******************')

simulation_name = Path(__file__).stem
print(f"Simulation name: {simulation_name}")


class wfc_controller:

    def __init__(self,
        n_turbines: int = 3,
        update_rate: float = 12.0,
        memory_size: int = 15,
        FAST_path_out: str = "/mnt/d/zero_model/fastfarm_simulation",
        min_yaw: float = -20.0,
        max_yaw: float = 20.0,
        base_prefix: str = simulation_name,
        input_path: str  = '/mnt/d/zero_model/floridyn_input_15013.json',
        output_path: str = '/mnt/d/zero_model/floridyn_output_15013.json',
        flag_path: str   = '/mnt/d/zero_model/floridyn_input_15013.ready'):

        self.n_turbines      = n_turbines
        self.t_update_rate   = update_rate
        self.base_prefix     = base_prefix
        self.min_yaw         = float(min_yaw)
        self.max_yaw         = float(max_yaw)
        self.input_path      = Path(input_path)
        self.output_path     = Path(output_path)
        self.flag_path       = Path(flag_path)        # NEW: flag file
        self._last_floridyn_t = -1e9

        self._t0               = None
        self._est_tick         = None
        self._last_est         = None
        self._last_logged_tick = None
        self._csv_idx          = 0

        self.base = Path(FAST_path_out)

        def _file_for_turbine(tid: int) -> Path:
            return self.base / f"{self.base_prefix}.T{tid}.out"
        self._file_for_turbine = _file_for_turbine

        self.history = [[dict() for _ in range(memory_size)] for _ in range(n_turbines)]
        self.opt_yaw = np.zeros(n_turbines, dtype=float)

        print('[FLORIDyn] Controller ready — waiting for MATLAB on Windows.')

    # ---- helpers -------------------------------------------------------------

    def _push_history(self, tidx: int, meas: dict):
        buf = self.history[tidx]
        buf[:-1] = buf[1:]
        buf[-1]  = meas

    def _call_floridyn(self, WS, WD, Thrust, Power, TI):
        """
        Send signals to MATLAB FLORIDyn and wait for gamma_opt.
        Uses flag-file protocol to avoid race conditions:
          1. Remove previous output
          2. Write JSON input
          3. Create .ready flag  <-- MATLAB waits for this
          4. Wait for output JSON from MATLAB
        """
        # 1. Remove stale output
        if self.output_path.exists():
            self.output_path.unlink()
        if self.flag_path.exists():
            self.flag_path.unlink()

        # 2. Write input JSON
        input_data = {
            'WS':     WS,           # scalar [m/s]
            'WD':     WD,           # scalar [deg, FLORIDyn convention]
            'Thrust': Thrust,       # [n_turb] [kN]
            'Power':  Power,        # [n_turb] [kW]
            'TI':     TI,           # [n_turb] [-]
        }
        with open(self.input_path, 'w') as f:
            json.dump(input_data, f)

        # 3. Create flag AFTER JSON is fully written
        self.flag_path.touch()

        # 4. Wait for MATLAB output (max 120s)
        for _ in range(120000):
            if self.output_path.exists():
                time_mod.sleep(0.05)
                try:
                    with open(self.output_path, 'r') as f:
                        output = json.load(f)
                    self.output_path.unlink()   # clean up
                    return output['gamma_opt']
                except (json.JSONDecodeError, KeyError):
                    time_mod.sleep(0.1)
                    continue
            time_mod.sleep(0.01)

        raise TimeoutError('[FLORIDyn] No response from MATLAB within 120s')

    def _average_FF_signals(self):
        """
        Read FAST.Farm .out files and return averaged signals.
        TI computed as mean(std(Wind_iX) / mean(Wind_iX)) over 9 wind points per turbine.
        """
        n_turb = self.n_turbines
        fs     = 10.0
        cutoff = 0.02
        order  = 2
        N      = 500    # last N samples (~50s at 10 Hz)
        N_eff  = N

        # --- Read header from T1 ---
        header_source = self._file_for_turbine(1)
        if not header_source.exists():
            raise FileNotFoundError(f"T1.out not found in {self.base}")

        with open(header_source, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        header_idx = None
        if len(lines) > 6 and ("Time" in lines[6] or "GenPwr" in lines[6]):
            header_idx = 6
        if header_idx is None:
            for i, ln in enumerate(lines[:200]):
                if ln.strip().startswith("Time"):
                    header_idx = i; break
        if header_idx is None:
            for i, ln in enumerate(lines[:200]):
                if "GenPwr" in ln:
                    header_idx = i; break
        if header_idx is None:
            raise ValueError("Cannot find header row")

        raw_hdr     = lines[header_idx].strip()
        header_cols = raw_hdr.split("\t") if "\t" in raw_hdr else re.split(r"\s+", raw_hdr)

        units_line   = lines[header_idx+1].strip() if header_idx+1 < len(lines) else ""
        looks_units  = any(tok in units_line for tok in ("(s)", "(m/s)", "(", ")", "[", "]", "\t"))
        data_start   = header_idx + 2 if looks_units else header_idx + 1

        def find_idx(cands, cols):
            for c in cands:
                if c in cols:
                    return cols.index(c)
            return None

        idx_genpwr    = find_idx(["GenPwr"],       header_cols)
        idx_genthrust = find_idx(["RotThrust"],    header_cols)
        idx_Ct        = find_idx(["RtFldCt"],      header_cols)
        idx_yaw       = find_idx(["YawPzn"],       header_cols)
        idx_windhub   = find_idx(["WindHubAngXY"], header_cols)
        idx_diskX     = find_idx(["WindDiskVelX"], header_cols)
        idx_diskY     = find_idx(["WindDiskVelY"], header_cols)
        idx_diskZ     = find_idx(["WindDiskVelZ"], header_cols)
        idx_Vel1X     = find_idx(["Wind1VelX"],    header_cols)
        idx_Vel2X     = find_idx(["Wind2VelX"],    header_cols)
        idx_Vel3X     = find_idx(["Wind3VelX"],    header_cols)
        idx_Vel4X     = find_idx(["Wind4VelX"],    header_cols)
        idx_Vel5X     = find_idx(["Wind5VelX"],    header_cols)
        idx_Vel6X     = find_idx(["Wind6VelX"],    header_cols)
        idx_Vel7X     = find_idx(["Wind7VelX"],    header_cols)
        idx_Vel8X     = find_idx(["Wind8VelX"],    header_cols)
        idx_Vel9X     = find_idx(["Wind9VelX"],    header_cols)

        if idx_genpwr is None or idx_windhub is None:
            raise ValueError(f"Essential columns missing: {header_cols}")

        def read_lastN(fp: Path, col_idx, n: int) -> np.ndarray:
            if col_idx is None:
                return np.full(n, np.nan)
            df = pd.read_csv(fp, sep="\t", skiprows=data_start, header=None, engine="c")
            if df.shape[1] <= col_idx:
                df = pd.read_csv(fp, delim_whitespace=True, skiprows=data_start,
                                 header=None, engine="c")
            if df.shape[1] <= col_idx:
                return np.array([])
            s = pd.to_numeric(df.iloc[:, col_idx], errors="coerce").dropna().to_numpy()
            return s[-n:] if s.size else np.array([])

        # --- Collect per-turbine signals ---
        P_list, T_list, Ct_list, yaw_list, WD_list = [], [], [], [], []
        diskX_list, diskY_list, diskZ_list = [], [], []
        Vel1X_list, Vel2X_list, Vel3X_list = [], [], []
        Vel4X_list, Vel5X_list, Vel6X_list = [], [], []
        Vel7X_list, Vel8X_list, Vel9X_list = [], [], []

        for tid in range(1, n_turb + 1):
            fp = self._file_for_turbine(tid)
            if not fp.exists():
                raise FileNotFoundError(f"Missing file: {fp.name}")
            P_list.append(read_lastN(fp, idx_genpwr,    N_eff))
            T_list.append(read_lastN(fp, idx_genthrust, N_eff))
            Ct_list.append(read_lastN(fp, idx_Ct,       N_eff))
            yaw_list.append(read_lastN(fp, idx_yaw,     N_eff))
            WD_list.append(read_lastN(fp, idx_windhub,  N_eff))
            diskX_list.append(read_lastN(fp, idx_diskX, N_eff))
            diskY_list.append(read_lastN(fp, idx_diskY, N_eff))
            diskZ_list.append(read_lastN(fp, idx_diskZ, N_eff))
            Vel1X_list.append(read_lastN(fp, idx_Vel1X, N_eff))
            Vel2X_list.append(read_lastN(fp, idx_Vel2X, N_eff))
            Vel3X_list.append(read_lastN(fp, idx_Vel3X, N_eff))
            Vel4X_list.append(read_lastN(fp, idx_Vel4X, N_eff))
            Vel5X_list.append(read_lastN(fp, idx_Vel5X, N_eff))
            Vel6X_list.append(read_lastN(fp, idx_Vel6X, N_eff))
            Vel7X_list.append(read_lastN(fp, idx_Vel7X, N_eff))
            Vel8X_list.append(read_lastN(fp, idx_Vel8X, N_eff))
            Vel9X_list.append(read_lastN(fp, idx_Vel9X, N_eff))

        def stack_matrix(lst):
            cols = []
            for s in lst:
                s = np.asarray(s, float).reshape(-1)
                if s.size >= N_eff:
                    cols.append(s[-N_eff:])
                elif s.size == 0:
                    cols.append(np.full(N_eff, np.nan))
                else:
                    cols.append(np.pad(s, (N_eff - s.size, 0), mode="edge"))
            return np.column_stack(cols)   # [N_eff x n_turb]

        GenPwr_raw  = stack_matrix(P_list)
        Thrust_raw  = stack_matrix(T_list)
        Ct_raw      = stack_matrix(Ct_list)
        Yaw_raw     = stack_matrix(yaw_list)
        WindDir_raw = stack_matrix(WD_list)
        diskX_raw   = stack_matrix(diskX_list)
        diskY_raw   = stack_matrix(diskY_list)
        diskZ_raw   = stack_matrix(diskZ_list)
        Vel1X_raw   = stack_matrix(Vel1X_list)
        Vel2X_raw   = stack_matrix(Vel2X_list)
        Vel3X_raw   = stack_matrix(Vel3X_list)
        Vel4X_raw   = stack_matrix(Vel4X_list)
        Vel5X_raw   = stack_matrix(Vel5X_list)
        Vel6X_raw   = stack_matrix(Vel6X_list)
        Vel7X_raw   = stack_matrix(Vel7X_list)
        Vel8X_raw   = stack_matrix(Vel8X_list)
        Vel9X_raw   = stack_matrix(Vel9X_list)

        # --- Low-pass filter (NOT applied to wind probes: filter removes turbulence) ---
        b, a   = butter(order, cutoff / (0.5 * fs), btype='low')
        padlen = 3 * (max(len(a), len(b)) - 1)
        def filt(M):
            return filtfilt(b, a, M, axis=0) if M.shape[0] > padlen else M

        GenPwr  = filt(GenPwr_raw)
        Thrust  = filt(Thrust_raw)
        Ct_f    = filt(Ct_raw)
        Yaw_f   = filt(Yaw_raw)
        WindDir = filt(WindDir_raw)
        diskX_f = filt(diskX_raw)
        diskY_f = filt(diskY_raw)
        diskZ_f = filt(diskZ_raw)

        sl = slice(-N_eff, None)

        # --- Wind speed (disk resultant, RMS) ---
        WindTot_vec = []
        for i in range(n_turb):
            disk_res = np.sqrt(diskX_f[:, i]**2 + diskY_f[:, i]**2 + diskZ_f[:, i]**2)
            WindTot_vec.append(float(np.sqrt(np.nanmean(disk_res[sl]**2))))
        WindTot_scalar = WindTot_vec[0]

        # --- Wind direction (T1 hub angle, mean) ---
        WindDir_scalar = float(np.nanmean(WindDir[sl, 0]))

        # --- Ct and Yaw averages ---
        Ct_vec  = np.nanmean(Ct_f[sl, :],  axis=0).tolist()
        Yaw_vec = np.nanmean(Yaw_f[sl, :], axis=0).tolist()

        # --- Thrust from Ct ---
        Thrust_vec = np.array([
            0.5 * 1.225 * ws**2 * (np.pi * 240.0**2) / 4
            * np.cos(np.deg2rad(WindDir_scalar - yaw))**1.88 * Ct
            for ws, yaw, Ct in zip(WindTot_vec, Yaw_vec, Ct_vec)
        ]) / 1e3
        Thrust_vec = Thrust_vec.tolist()

        # --- Power [kW] ---
        Power_vec = np.nanmean(GenPwr[sl, :], axis=0).tolist()

        # --- TI per turbine from 9 wind probe X-components (RAW, not filtered) ---
        # TI_k = std(VelkX) / mean(VelkX)  for each probe k=1..9
        # TI_turb = cubic mean of TI_1..TI_9
        all_Vel_raw = [Vel1X_raw, Vel2X_raw, Vel3X_raw,
                       Vel4X_raw, Vel5X_raw, Vel6X_raw,
                       Vel7X_raw, Vel8X_raw, Vel9X_raw]

        TI_vec = []
        for i in range(n_turb):
            ti_probes = []
            for Velk_raw in all_Vel_raw:
                probe = Velk_raw[sl, i]
                m = np.nanmean(probe)
                s = np.nanstd(probe)
                if m > 0.1:
                    ti_probes.append(s / m)
            if ti_probes:
                # cubic mean: (mean(TI^3))^(1/3)
                TI_i = float(np.cbrt(np.nanmean(np.array(ti_probes)**3)))
            else:
                TI_i = 0.06
            TI_vec.append(TI_i)

        # --- Console print ---
        print(f'  WindDir (T1): {WindDir_scalar:.2f} deg')
        for i in range(n_turb):
            print(f'  T{i+1}: WS={WindTot_vec[i]:.2f} m/s  '
                  f'Thrust={Thrust_vec[i]:.1f} kN  '
                  f'Power={Power_vec[i]:.0f} kW  '
                  f'TI={TI_vec[i]*100:.2f}%  '
                  f'Yaw={Yaw_vec[i]:.2f} deg  '
                  f'Ct={Ct_vec[i]:.4f}')

        return WindDir_scalar, Thrust_vec, WindTot_scalar, Power_vec, TI_vec

    def _ensure_estimate(self, current_time: float):
        if self._t0 is None:
            self._t0 = float(current_time)

        tick = int(round(current_time - self._t0))

        if self._est_tick == tick and self._last_est is not None:
            return self._last_est

        signals       = self._average_FF_signals()
        self._last_est = signals
        self._est_tick = tick
        return signals

    # ---- ROSCO interface -----------------------------------------------------

    def update_setpoints(self, id: int, current_time: float, measurements: dict):
        print(f'[DEBUG] id={id}, t={current_time}')
        tidx = int(id) - 1
        if tidx < 0 or tidx >= self.n_turbines:
            return {"ZMQ_YawOffset": 0.0, "ZMQ_PitOffset(1)": 0.0,
                    "ZMQ_PitOffset(2)": 0.0, "ZMQ_PitOffset(3)": 0.0}

        self._push_history(tidx, measurements or {})

        if (id == 1 and float(current_time) >= 60.0
                and (current_time - self._last_floridyn_t) >= self.t_update_rate):
            self._last_floridyn_t = current_time
            signals = self._ensure_estimate(float(current_time))
            
            if signals is not None:
                WindDir_scalar, Thrust_vec, WindTot_scalar, Power_vec, TI_vec = signals

                # Conversione prima del print
                WindDir_floridyn = -WindDir_scalar + 270.0

                print(f'[FLORIDyn INPUT] t={current_time:.1f}s')
                print(f'  WS={WindTot_scalar:.3f} m/s  WD_FF={WindDir_scalar:.2f} deg  WD_FLORIDyn={WindDir_floridyn:.2f} deg')
                for i in range(self.n_turbines):
                    print(f'  T{i+1}: Power={Power_vec[i]:.1f} kW  Thrust={Thrust_vec[i]:.1f} kN  TI={TI_vec[i]*100:.2f}%')

                gamma_opt = self._call_floridyn(
                    WindTot_scalar,
                    WindDir_floridyn,
                    Thrust_vec,
                    Power_vec,
                    TI_vec
                )
                yaws = np.array(gamma_opt).ravel()
                self.opt_yaw = np.clip(yaws, self.min_yaw, self.max_yaw)
                print(f'[FLORIDyn] yaw setpoints: {np.round(self.opt_yaw, 2).tolist()}')

        yaw_FLORIDyn = float(self.opt_yaw[tidx]) if tidx < len(self.opt_yaw) else 0.0

        return {"ZMQ_YawOffset":    yaw_FLORIDyn,
                "ZMQ_PitOffset(1)": 0.0,
                "ZMQ_PitOffset(2)": 0.0,
                "ZMQ_PitOffset(3)": 0.0}


if __name__ == "__main__":
    here    = os.path.dirname(os.path.abspath(__file__))
    logfile = os.path.join(here, "wfc_FLORIDyn.log")
    addr    = "tcp://localhost:15013"

    server = wfc_zmq_server(addr, timeout=5000.0, verbose=True, logfile=logfile)

    ctrl = wfc_controller(
        n_turbines=3,
        update_rate=12.0,
        memory_size=15,
        FAST_path_out="/mnt/d/zero_model/fastfarm_simulation",
        min_yaw=-20.0,
        max_yaw=20.0,
        base_prefix="15013",
        input_path='/mnt/d/zero_model/floridyn_input_15013.json',
        output_path='/mnt/d/zero_model/floridyn_output_15013.json',
        flag_path='/mnt/d/zero_model/floridyn_input_15013.ready',
    )
    server.wfc_controller = ctrl
    server.runserver()