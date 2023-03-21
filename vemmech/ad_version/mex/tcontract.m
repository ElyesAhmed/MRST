function varargout = tcontract(varargin)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

   filename = 'tcontract.cpp';
   INCLUDE = {};
   OPTS= {'-O'};
   SRC = {filename};
   [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags();
   
   %@@ DEBUG
   %CXXFLAGS = {'CXXFLAGS=$CXXFLAGS -D_GNU_SOURCE -D_GLIBCXX_DEBUG=1 -g -fPIC -O3 -std=c++11 -ffast-math -march=native -fopenmp'}
   
   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});
   
   [varargout{1:nargout}] = tcontract(varargin{:});
end
