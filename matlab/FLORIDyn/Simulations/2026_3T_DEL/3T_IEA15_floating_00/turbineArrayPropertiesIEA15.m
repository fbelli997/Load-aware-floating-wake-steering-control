% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

function tProp = turbineArrayPropertiesIEA15()
%% Return location and type of the turbines

% % Foot position of the turbine x,y,z (NOT nacelle position, Ground level!)
% tProp.Pos = [...
%             300  1860 0;...     % T0
%             1500 1860 0;...     % T1
%             2700 1860 0;...     % T2
%             ]; 

% % Wind direction fixed on FASTfarm
tProp.Pos = [...
            300  1650 0;...     % T0
            1500 1790 0;...     % T1
            2700 1930 0;...     % T2
            ]; 

% tProp.Pos = [...
%             300  1930 0;...     % T0
%             1500 1790 0;...     % T1
%             2700 1650 0;...     % T2
%             ]; 


% Turbine type, the parameters are stored centralized in
% Data/TurbineData/getTurbineData.m
tProp.Type = {...
            'IEA 15MW',...
            'IEA 15MW',...
            'IEA 15MW'};


% Initial world orientation of the Wind Turbine. Direction is 
% Orthogonal for wind from
% Upload Cp Ct LUT from FASTfarm and predictors and correctors for loads

% %%
% % ----------------------------
% % 1) Python environment
% % ----------------------------
% pe = pyenv;
% if pe.Status == "NotLoaded"
%     pyenv("Version", "C:\Users\fbell\AppData\Local\Programs\Python\Python311\python.exe");
% end
% disp(pyenv)
% 
% % Import Python libraries
% joblib_mod = py.importlib.import_module("joblib");
% torch_mod  = py.importlib.import_module("torch");
% 
% % ----------------------------
% % 2) Paths
% % ----------------------------
% joblib_path = "gpr_residual_tower.joblib";
% pt_path     = "del_tower_uncertainty_grouped.pt";
tProp.VCpCt         = load("VCpCt_15MW_FF.mat");
tProp.MatFloat      = load("Activefloat.mat");
tProp.MatStateSpace = load("D:\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00\LinModel_3DOF_Augmented_Wave.mat");
tProp.ptfm_initial  = deg2rad(-1.5);

% load corrector and predictor:
% tProp.GPR_corrector_tower = load("gpr_residual_tower.mat");
tProp.NN_predictor_tower  = load("del_tower_uncertainty_grouped_wind_dir.mat");
tProp.NN_corrector_tower = load("nn_residual_tower_wind_dir.mat");
% tProp.GPR_corrector_blade = load("gpr_residual_blade.mat");
tProp.NN_corrector_blade = load("nn_residual_blade_wind_dir.mat");
tProp.NN_predictor_blade  = load("del_blade_uncertainty_grouped_wind_dir.mat");

%% Has to fit the turbine states structure of the selected FLORIS model (see
% States.m), Gaussian: 'a','yaw'(deg,counter-clockwise,relative),'TI_added'
tProp.Init_States = [...
    0.33 17 0; ...
    0.33 17 0; ...
    0.33 0 0; ...
    ];


% TI here is hardcoded, you should think about having a relation with the
% TI calculated on the excel file:"D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00\Data\WindTIConstant.csv"

% NOTE ON THE TURBINE ORIENTATION, THE YAW ANGLE AND THE WIND DIRECTION
%   The wind direction and the turbine orientation are clockwise defined,
%   the yaw angle is the result of 
%       wind direction - orientation = yaw angle
%   Thus, the yaw angle is defined in clockwise direction.
%   Wind from -x to +x has the wind direction 270 deg.
%   Wind from +y to -y has the wind direction 0 / 360 deg.
end