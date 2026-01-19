classdef MolecDiffusFlowDiscretization < FlowDiscretization
    %Discretization and state function grouping for bio-chemistry flow
    % Author: [Stéphanie Delage Santacreu]
    % Date: [16/09/2025]
    % Organization: [Université de Pau et des Pays de l'Adour, E2S UPPA, CNRS, LFCR, UMR5150, Pau, France]
    % ---------------------------------------------------------------------------


    properties
        MolecularDiffPhaseFlux
    end

    methods
        %-----------------------------------------------------------------%
        function props = MolecDiffusFlowDiscretization(model)
            % Inherit most of the state functions from FLuxDiscretization
            props = props@FlowDiscretization(model);
            % Fluid flow transmissibility
            if model.dynamicFlowTrans
                % Dynamic transmissibility. Conductivities are caluclated
                % at each nonlinear iteration, and used to compute
                % corresponding transmissibilities
                props = props.setStateFunction('Transmissibility', ...
                    DynamicFlowTransmissibility(model, 'Permeability'));
            else
                % Static transmissibilities already been set up by parent
            end

            if model.dynamicFlowPv
                % Dynamic transmissibility. Conductivities are caluclated
                % at each nonlinear iteration, and used to compute
                % corresponding transmissibilities
                props = props.setStateFunction('PoreVolume', ...
                    DynamicFlowPoreVolume(model, 'Porosity'));
            else
                % Static transmissibilities already been set up by parent
            end
            %molecular diffusion
            props = props.setStateFunction('MolecularDiffPhaseFlux', ComponentMolecularDiffPhaseFlux(model));


        end

        %-----------------------------------------------------------------%
        function [acc, flux, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            [acc, flux, names, types] = componentConservationEquations@FlowDiscretization(fd, model, state, state0, dt);
            flowState = fd.buildFlowState(model, state, state0, dt);
            act = model.getActivePhases();
            ncomp = model.getNumberOfComponents();
            nph = sum(act);
            if model.moleculardiffusion
                Jmoldiff = model.getProps(flowState, 'MolecularDiffPhaseFlux');
                for c = 1:ncomp
                    for ph = 1:nph
                        flux{c} = flux{c} + Jmoldiff{c,ph};
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