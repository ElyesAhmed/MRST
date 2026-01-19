classdef ComponentPhaseMolecularDiffFlux < StateFunction
    % Flux of each component, in each phase
    % Molecular diffusion flux : 
    % Quirk-Millington model for both phases to evaluate the diffusion
    % coefficients.
    % For the gas phase, the binary diffusion coefficients are estimated 
    % thanks to the  Stefan-Maxwell equation simplified by Blanc's law. 
    % For the liquid phase, the binary diffusion coefficients comes from
    % the literature.
    % Author: [Stéphanie Delage Santacreu]
    % Date: [16/09/2025]
    % Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
    % ---------------------------------------------------------------------------

    properties
        
    end

    methods
        function gp = ComponentPhaseMolecularDiffFlux(model, varargin)
            gp@StateFunction(model);                
            gp = gp.dependsOn('Density', 'PVTPropertyFunctions');
            gp = gp.dependsOn({'PoreVolume'}, 'PVTPropertyFunctions');%sds modif
             gp = gp.dependsOn('s', 'state');
             gp = gp.dependsOn('x', 'state');
             gp = gp.dependsOn('y', 'state');
             gp = gp.dependsOn('pressure', 'state');
             gp = gp.dependsOn('T', 'state');

            gp.label = 'J_{i,\alpha}';
        end

        function J = evaluateOnDomain(prop, model, state)            
            ncomp = model.getNumberOfComponents;            
            nph = model.getNumberOfPhases;
            J = cell(ncomp, nph);
            J = cellfun(@(x) 0, J, 'UniformOutput', false);      
            if isfield(state,'x')
                nm = model.getPhaseNames();
                namecp = model.compFluid.names();
                indices = struct('H2', find(strcmp(namecp, 'H2')), ...
                    'C1', find(strcmp(namecp, 'C1')), ...
                    'CO2', find(strcmp(namecp, 'CO2')), ...
                    'H2O', find(strcmp(namecp, 'H2O')), ...
                    'N2', find(strcmp(namecp, 'N2')), ...
                    'C2', find(strcmp(namecp, 'C2')), ...
                    'C3', find(strcmp(namecp, 'C3')), ...
                    'H2S', find(strcmp(namecp, 'H2S')), ...
                    'NC4', find(strcmp(namecp, 'NC4')));

                fields = fieldnames(indices);
                nfields=numel(fields);
                mol_diff=zeros(nfields,2);
                param_LJ=zeros(nfields,2);

                %At  25 degrees Celsius: coef dans l'eau douce (coef not used in gas)
                %Attention: salinity reduces coef from 10 to 30 %
                % coeffs = struct(...
                %     'H2',  [4.5e-9, 6.1e-5], ...
                %     'C1', [1.5e-9, 1.6e-5], ...
                %     'H2O', [2.3e-9, 1.5e-5], ...
                %     'CO2', [1.9e-9, 1.4e-5], ...
                %     'N2',  [2.0e-9, 1.8e-5], ...
                %     'C2',  [1.2e-9, 2.5e-5], ...
                %     'C3',  [1.0e-9, 2.2e-5], ...
                %     'H2S',  [1.5e-9, 2.2e-5], ...
                %     'NC4', [0.8e-9, 1.9e-5]);

                %At  40 degrees Celsius: coef dans l'eau douce (coef not used in gas)
                %Attention: salinity reduces coef from 10 to 30 %
                coeffs = struct(...
                    'H2',  [6.44e-9, 6.1e-5], ...
                    'C1', [2.15e-9, 1.6e-5], ...
                    'H2O', [3.29e-9, 1.5e-5], ...
                    'CO2', [2.72e-9, 1.4e-5], ...
                    'N2',  [2.86e-9, 1.8e-5], ...
                    'C2',  [1.72e-9, 2.5e-5], ...
                    'C3',  [1.43e-9, 2.2e-5], ...
                    'H2S',  [2.15e-9, 2.2e-5], ...
                    'NC4', [1.15e-9, 1.9e-5]);


                coeffs_LJ=struct(... %paramètres de Lennard-Jones (diametre moyen (Angstrom), potentiel(K))
                    'H2',  [2.92,59.7], ...
                    'C1',  [3.758,148.6], ...
                    'H2O',  [2.641,809.1], ...
                    'CO2',  [3.996,195.2], ...
                    'N2',   [3.798,71.4], ...
                    'C2',   [4.443,215.7], ...
                    'C3',   [5.118,237.1], ...
                    'H2S',   [3.60,301.0], ...
                    'NC4', [5.206,289.5]);

                for i = 1:nfields
                    comp = fields{i};
                    indComp = indices.(comp);

                    if ~isempty(indComp) && isfield(coeffs, comp)
                        mol_diff(indices.(comp),:)= coeffs.(comp);
                        param_LJ(indices.(comp),:)=coeffs_LJ.(comp);
                    end
                end

               rho = prop.getEvaluatedExternals(model, state, 'Density'); 
               [p, T] = model.getProps(state, 'pressure', 'temperature'); 
               coeff1=T.^(1.5)./(9.869e-6.*p); %1.e5.*p en atm et p en Pa 
               Molmass=1.e3.*model.compFluid.molarMass;
               SigLJ=param_LJ(:,1);
               EpsLJ=param_LJ(:,2);

               avg = model.operators.faceAvg;

               %bioclogging
               if model.dynamicFlowPv()
                   if model.bacteriamodel
                       poro=model.rock.poro(p,state.nbact);
                   else
                       poro=model.rock.poro(p);
                   end
               else
                   poro= model.rock.poro;
               end


               
               L_ix = model.getLiquidIndex();
               V_ix = model.getVaporIndex();
                
               % Define diffusion coefficients in m²/s for liquid and gas phases
               sqrtMij=zeros(ncomp,ncomp);
               sqrtEpsij=zeros(ncomp,ncomp);
               Sigij2=zeros(ncomp,ncomp);
               for c = 1:ncomp
                   for cj = 1:ncomp
                       sqrtMij(c,cj)=sqrt(2*Molmass(c)*Molmass(cj)/(Molmass(c)+Molmass(cj)));
                       sqrtEpsij(c,cj)=sqrt(EpsLJ(c)*EpsLJ(cj));
                       Sigij2(c,cj)=0.25*(SigLJ(c)+SigLJ(cj))*(SigLJ(c)+SigLJ(cj));
                   end
               end
               Dij=1.e-4.*0.001858./(sqrtMij.*Sigij2); %en m2/s


               for c = 1:ncomp
                   if iscell(state.x)
                       xc = state.x{c};
                   else
                       xc = state.x(c);
                   end
                   if iscell(state.y)
                       yc = state.y{c};
                   else
                       yc = state.y(c);
                   end
                   for ph = 1:nph
                       s = model.getProp(state, ['s', nm(ph)]);
                       tau_mq=(s.*poro).^(7/3).*poro.^(-2);
                       
                       if (ph==L_ix) 
                           D_diffl = avg(s.*rho{ph}.*mol_diff(c,ph).*tau_mq.*poro);%Millington and Quirk model
                           J{c, ph} = - D_diffl.*model.operators.Grad(xc);
                       elseif (ph==V_ix)
                           %calcul des Dij
                           D_diffij_inv =yc.*0;

                           for cj = 1:ncomp
                               if iscell(state.y)
                                   ycj = state.y{cj};
                               else
                                   ycj = state.y(cj);
                               end
                               if (cj~=c)
                                   D_diffij_inv =  D_diffij_inv  +ycj./Dij(c,cj);
                               end
                           end
                           D_diffij=coeff1./max(D_diffij_inv,1.e-16);
                           D_diff = avg(s.*rho{ph}.*D_diffij.*tau_mq.*poro);%Millington and Quirk model
                           J{c, ph} = - D_diff.*model.operators.Grad(yc);
                       end
                   end
               end
            end
        end
    end
end
%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}