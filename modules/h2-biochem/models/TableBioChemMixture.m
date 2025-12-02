classdef TableBioChemMixture
    % Create MRST parameters lis for biochemical reaction from name
    properties
        metabolicReaction % name of bioreactions :'MethanogenicArchae','SulfateReducingBacteria'
        rH2 %  Names of bioreactions reactants
        rsub
        pH2O%  Names of bioreactions products
        p2
        gamrH2 %  Associated Stoichiometric coefficients
        gamrsub
        gampH2O
        gamp2
        Y_H2 %  Yield coefficient
        alphaH2 %  Half-saturation constant for H2
        alphasub %  Half-saturation constant for substrate (CO2 or SO4)
        Psigrowthmax %  Maximum growth rate
        b_bact %  Constant decay rate
        nbactMax  %  Maximum number of microorganisms per volume
    end

    methods
        function biochemmix = TableBioChemMixture(metabolicReactions)
            % Create fluid biochemical reaction from metabolicReaction name
            %
            % SYNOPSIS:
            %   f = TableBioChemMixture({'MethanogenicArchae'});
            %
            % PARAMETERS:
            %   names    - Cell array of valid names. See `getFluidList` for
            %              valid names.
            %              (optional).
            % RETURNS:
            %   biochemmix - Initialized compositional mixture.
            if ischar(metabolicReactions)
                metabolicReactions = {metabolicReactions};
            end

            nbioreact = numel(metabolicReactions);
            [metabolicReaction,rH2,rsub,pH2O,p2]=deal(strings(1,nbioreact));
            [gamrH2,gamrsub,gampH2O,gamp2,Y_H2,alphaH2,alphasub,...
                Psigrowthmax,b_bact,nbactMax] = deal(zeros(1, nbioreact));


            biochemfluids = bioChemFluidsStructs();
            validChoices = TableBioChemMixture.getFluidList();

            ok = ismember(lower(metabolicReactions), lower(validChoices));

            if ~all(ok)
                s ='Unable to create bioreaction paramaters. The following names were not known: ';
                msg = [s, sprintf('%s ', metabolicReactions{~ok})];
                error(msg);
            end

            for i = 1:numel(metabolicReactions)
                isF = strcmpi(metabolicReactions{i}, validChoices);
                str = biochemfluids(isF);
                metabolicReaction(i)=str.metabolicReaction;
                rH2(i)=str.rH2;
                rsub(i)=str.rsub;
                pH2O(i)=str.pH2O;
                p2(i)=str.p2;
                gamrH2(i) = str.gamrH2;
                gamrsub(i) = str.gamrsub;
                gampH2O(i) = str.gampH2O;
                gamp2(i) = str.gamp2;
                Y_H2(i) = str.Y_H2;
                alphaH2(i) = str.alphaH2;
                alphasub(i) = str.alphasub;
                Psigrowthmax(i) = str.Psigrowthmax;
                b_bact(i) = str.b_bact;
                nbactMax(i) = str.nbactMax;
            end



            biochemmix.metabolicReaction=metabolicReaction;
            biochemmix.rH2=rH2 ;
            biochemmix.rsub=rsub;
            biochemmix.pH2O=pH2O;
            biochemmix.p2=p2;
            biochemmix.gamrH2=gamrH2;
            biochemmix.gamrsub=gamrsub;
            biochemmix.gampH2O=gampH2O;
            biochemmix.gamp2=gamp2;
            biochemmix.Y_H2=Y_H2;
            biochemmix.alphaH2=alphaH2;
            biochemmix.alphasub=alphasub ;
            biochemmix.Psigrowthmax=Psigrowthmax;
            biochemmix.b_bact=b_bact ;
            biochemmix.nbactMax=nbactMax;




        end
    end
    methods (Static)
        function varargout = getFluidList()
            biochemfluids = bioChemFluidsStructs();
            metabolicReactions = {biochemfluids.metabolicReaction};
            if nargout > 0
                varargout{1} = metabolicReactions;
            else
                disp('Possible biochemfluids choices are:');
                fprintf('%s\n', metabolicReactions{:});
            end
        end
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
