classdef GrowthBactRateSRC < StateFunction
    % Bacterial growth rate computation for compositional simulations
    %
    % SYNOPSIS:
    %   gr = GrowthBactRateSRC(model, 'property1', value1, ...)
    %
    % DESCRIPTION:
    %   Computes the specific growth rate coefficient (kinetic term only)
    %   for each microbial population using Monod kinetics.
    %
    %   For methanogens and acetogens, H2 and CO2 are taken from the
    %   liquid-phase EOS composition (x). For sulfate reducers (SRB),
    %   the substrate SO4 is taken from the aqueous tracer `tracerSO4`
    %   and converted to mole fraction using the liquid molar density.

    properties
        % No additional properties
    end

    methods
        function gp = GrowthBactRateSRC(model, varargin)
            % Constructor
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('x', 'state');
            % If SRB is active, we need the tracer and liquid Z-factor
            if isprop(model, 'sulfateReduction') && model.sulfateReduction
                gp = gp.dependsOn('tracerSO4', 'state');
                gp = gp.dependsOn('Z_L', 'state');
            end
            gp.label = '\Psi_{growth}';
        end

        function Psigrowth = evaluateOnDomain(prop, model, state)
            % Compute specific growth rate coefficient [1/s]
            %
            % Returns: Psigrowthmax * axH2 * axsub [1/s]
            % Used with BacterialMass: source = Psigrowth * BacterialMass

            if isprop(model, 'ReservoirModel') && ~isempty(model.ReservoirModel)
                rm = model.ReservoirModel;
            else
                rm = model;
            end
            bcrm = rm.biochemFluid;
            namecp = rm.getComponentNames();
            nbioreact = bcrm.nbioreact;

            % Initialize output
            Psigrowth = cell(1, nbioreact);
            [Psigrowth{:}] = deal(0);

            % Get liquid mole fractions (EOS components)
            x = rm.getProp(state, 'x');

            % Loop over reactions
            for i = 1:nbioreact
                % Find H2 index (required for all reactions)
                idx_H2 = find(strcmpi(namecp, bcrm.rH2(i)), 1);
                if isempty(idx_H2)
                    continue;   % H2 not found – skip this reaction
                end
                % H2 mole fraction
                if iscell(x)
                    xH2 = x{idx_H2};
                else
                    xH2 = x(:, idx_H2);
                end

                % Handle substrate depending on reaction type
                if strcmp(bcrm.metabolicReaction(i), 'SulfateReducingBacteria')
                    % SRB: substrate is sulfate tracer (not in EOS)
                    if isfield(state, 'tracerSO4')
                        so4_conc = state.tracerSO4;   % mol/m3 liquid
                    else
                        so4_conc = zeros(size(xH2));
                    end
                    % Compute liquid molar density from Z_L
                    if isfield(state, 'Z_L')
                        Z_L = state.Z_L;
                    else
                        Z_L = ones(size(xH2));
                    end
                    R = 8.314;   % Pa·m3/(mol·K)
                    % pressure, T may be ADI objects; use value() if needed
                    P = state.pressure;
                    T = state.T;
                    rho_molar = P ./ (Z_L .* R .* T);   % mol/m3
                    % Convert tracer to mole fraction
                    xsub = so4_conc ./ rho_molar;
                    % Use the half‑saturation constant for sulfate (already in mol/mol)
                    alphasub = bcrm.alphasub(i);
                else
                    % Methanogens and acetogens: substrate is an EOS component
                    idx_sub = find(strcmpi(namecp, bcrm.rsub(i)), 1);
                    if isempty(idx_sub)
                        continue;
                    end
                    if iscell(x)
                        xsub = x{idx_sub};
                    else
                        xsub = x(:, idx_sub);
                    end
                    alphasub = bcrm.alphasub(i);
                end

                % Now compute Monod terms
                alphaH2 = bcrm.alphaH2(i);
                Psigrowthmax = bcrm.Psigrowthmax(i);

                axH2 = xH2 ./ (alphaH2 + xH2);
                axsub = xsub ./ (alphasub + xsub);

                Psigrowth{i} = Psigrowthmax .* axH2 .* axsub;
            end
        end
    end
end