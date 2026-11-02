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

function intOPs = interpolateOPs(nT, StartI, dep, States_OP, posBase, nOP)% T.nT, T.StartI, T.dep, T.States_OP, T.posBase, T.nOP
%INTERPOALTEOPS Interpolates the OPs at the closest location of the
%turbines which are then used to set up the temporary wind farms.

intOPs = cell(nT,1);
StartI = StartI;

for iT = 1:nT % Every turbine
    % [id OP1, rat OP1, id OP2, rat OP2]
    intOPs{iT} = zeros(length(dep{iT}),4);
    
    for iiT = 1:length(dep{iT}) % every influencing turbine
        % actual turbine index
        iiaT = dep{iT}(iiT);
        
        % distances of influencing OPs to turbine
        dist = sqrt(sum(...
            (States_OP(StartI(iiaT)+(0:nOP-1),1:2)...
            - posBase(iT,1:2)).^2,2));
        
        % get closest OP
        [~,i_min] = sort(dist);
        
        % get second closest
        if i_min(1) == 1
            % first OP (unlikely)
            intOPs{iT}(iiT,:) = [StartI(iiaT), 1, StartI(iiaT)+1, 0];
        elseif i_min(1) == nOP
            % last OP (possible)
            intOPs{iT}(iiT,:) = [StartI(iiaT)+nOP-2, 0, ...
                StartI(iiaT)+nOP-1, 1];
        else
            %% use two closest OPs to interpolate influence
            indOP1 = StartI(iiaT) - 1 + i_min(1);
            indOP2 = StartI(iiaT) - 1 + i_min(2);
            % Relative Distance: d = -((b-a)⊤(a-c))/((b-a)⊤(b-a))
            % Interpolation location: a+(b-a)d
            % a: OP close, b: OP far, c: turb_location
            % Interpolation: d⋅b + (1-d)a
            % Bound d to [0,1]
            a = States_OP(indOP1,1:2)';
            b = States_OP(indOP2,1:2)';
            c = posBase(iT,1:2)';
            
            d = ((b-a)'*(c-a))/((b-a)'*(b-a));
%             figure(1)
%             hold on
%             scatter(d,min(max(d,0),1));
%             hold off
%             pause(0.1)
            d = min(max(d,0),1);
            
            r1 = 1-d;
            r2 = d;
            
            % DEBUGGING FIGURE
%             figure
%             scatter(States_OP(indOP1,1),States_OP(indOP1,2))
%             hold on
%             scatter(States_OP(StartI(iiaT)+(0:nOP-1),1),...
%                 States_OP(StartI(iiaT)+(0:nOP-1),2))
%             scatter(States_OP(indOP2,1),States_OP(indOP2,2))
%             scatter(posBase(iT,1),posBase(iT,2))
%             scatter(States_OP(indOP1,1)*r1 + States_OP(indOP2,1)*r2,...
%                 States_OP(indOP1,2)*r1 + States_OP(indOP2,2)*r2,...
%                 'filled')
%             tpos = [posBase(iT,1),posBase(iT,2)];
%             viscircles([tpos;tpos;tpos],...
%                 [dist(i_min),dist(i_min-1),dist(i_min+1)]);
%             hold off
%             legend('OP1','OP2','Turb','InterpOP')
            intOPs{iT}(iiT,:) = [indOP1,r1,indOP2,r2];
        end
    end % every influencing turbine
end % every turbine

end

