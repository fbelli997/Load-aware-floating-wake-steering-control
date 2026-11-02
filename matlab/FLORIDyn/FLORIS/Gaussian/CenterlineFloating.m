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

function delta = CenterlineFloating(States_OP,States_T,States_WF,paramFLORIS,D,T, Ct)
%CENTERLINE model

%% Parameters
k_a     = paramFLORIS.k_a;
k_b     = paramFLORIS.k_b;
alpha   = paramFLORIS.alpha;
beta    = paramFLORIS.beta;

%% States
% C_T     = CalcCt(States_T(:,1),States_T(:,2)); % check the multiplication
% with yaw as FLORIS does. OLD METHOD FROM AXIAL INDUCTION factor.
nT  = numel(T);              
Ntot = size(States_OP, 1);   
Nop  = Ntot / nT;

assert(mod(Ntot, nT) == 0, 'Ntot must be divisible by nT');

C_T = repelem(Ct(:), Nop);  
yaw     = -deg2rad(States_T(:,2));
I       = sqrt(States_T(:,3).^2+States_WF(:,3).^2); % I_f & I_0
OPdw    = States_OP(:,4);
tilt    = repelem(T(:), Nop);
tilt_deg = rad2deg(tilt);

%% Implementation of true misalignment method
yaw_i = acos(cos(yaw).*cos(tilt));

%% Calc x_0 (Core length)
% [1] Eq. 7.3
x_0 = (cos(yaw_i).*(1+sqrt(1-C_T))./...
    (sqrt(2)*(alpha*I+beta*(1-sqrt(1-C_T))))).*D;

%% Calc k_z and k_y based on I
%[2] Eq.8
k_y = k_a*I + k_b;
k_z = k_y;

%% Get field width y
% [1] Eq. 7.2
% To fit the field width, the value linearly increases from 0 to max for dw
% positions before x_0
zs = zeros(size(OPdw));
sig_y = ...
    max([OPdw-x_0,zs],[],2)   .* k_y +...
    min([OPdw./x_0,zs+1],[],2).* cos(yaw_i) .* D/sqrt(8);

%% Get field width z
% [1] Eq. 7.2
sig_z = ...
    max([OPdw-x_0,zs],[],2)    .* k_z +...
    min([OPdw./x_0,zs+1],[],2) .* D/sqrt(8);

%% Calc Theta
%[1] Eq. 6.12
Theta = 0.3*yaw_i./cos(yaw_i).*(1-sqrt(1-C_T.*cos(yaw_i)));

%% Calc Delta / Deflection
%[1] Eq. 7.4 (multiplied with D and the second part disabled for dw<x_0)
%   Part 1 covering the linear near field and constant for the far field
delta_nfw  = Theta.*min([OPdw,x_0],[],2);

%   Part 2, smooth angle for the far field, disabled for the near field.
delta_fw_1 = Theta/14.7.*sqrt(cos(yaw_i)./(k_y.*k_z.*C_T)).*(2.9+1.3*sqrt(1-C_T)-C_T);
delta_fw_2 = log(...
    (1.6+sqrt(C_T)).*...
    (1.6.*sqrt((8*sig_y.*sig_z)./(D.^2.*cos(yaw_i)))-sqrt(C_T))./(...
    (1.6-sqrt(C_T)).*...
    (1.6.*sqrt((8*sig_y.*sig_z)./(D.^2.*cos(yaw_i)))+sqrt(C_T))...
    ));

% Combine deflection parts
delta_mis = delta_nfw + ...
    (sign(OPdw-x_0)/2+0.5).*...
    delta_fw_1.*delta_fw_2.*D;

% yaw, tilt, deflect_mis devono essere scalari (in gradi per yaw/tilt)
den = sqrt(sin(yaw).^2 + (sin(tilt).*cos(yaw)).^2);

y_PI = zeros(length(den), 2);   % 55x2, default 0
idx  = (den ~= 0);

y_PI(idx,:) = [ ...
    sin(yaw(idx)), ...
    sin(tilt(idx)).*cos(yaw(idx)) ...
] ./ den(idx);

deflect_yz = delta_mis .* y_PI;   


deltaZ = deflect_yz(:,2);
deltaY = -deflect_yz(:,1); 


% Deflection in y and z direction
delta = [deltaY, deltaZ];
end

