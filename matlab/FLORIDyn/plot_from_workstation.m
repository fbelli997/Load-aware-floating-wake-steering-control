% 1. Aggiungi i path delle funzioni (stessa cosa che fai sulla workstation)
addPaths;
addFLORISPaths;
addFLORIDynPaths;

load("C:\Users\fbell\OneDrive - Politecnico di Milano\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00\Results\results_CLC_RA.mat")
Sim.PathToSim = [pwd filesep];  % usa la cartella corrente del laptop
set(0, 'DefaultFigureVisible', 'on');
plot_CLC_results(time_hist, time_hist_control, ws_hist, wd_hist, ws_std_hist, wd_std_hist, ...
                 yaw_hist, DEL_NN_hist_blade, DEL_NN_hist_tower, ...
                 DD_hist, YD_hist, TI_hist, nds_hist, T, Sim);