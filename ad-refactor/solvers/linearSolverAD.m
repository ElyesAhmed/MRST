classdef linearSolverAD < handle
    % Base class for a nonlinear solver
   properties
       tolerance
       maxIterations
       extraReport
       
       verbose
   end
   methods
       function solver = linearSolverAD(varargin)
           solver.tolerance     = 1e-8;
           solver.maxIterations = 25;
           solver.extraReport   = false;
           solver.verbose       = mrstVerbose();
           if numel(varargin)
              solver = mergeSupportedOptions(solver, varargin{:});
           end
       end
       
       function [result, report] = solveLinearSystem(solver, A, b) %#ok
           % Solve the linear system to a given tolerance
           error('Superclass not meant for direct use')
       end
       
       function [dx, result, report] = solveLinearProblem(solver, problem, model)
           % Solve a linearized problem
           problem = problem.assembleSystem();
           
           timer = tic();
           [result, report] = solver.solveLinearSystem(problem.A, problem.b); 
           report.SolverTime = toc(timer);
           
           dx = solver.storeIncrements(problem, result);
       end
       
       function dx = storeIncrements(solver, problem, result) %#ok
           % Store the increments seperately - copy/paste from SolveEqsADI
           numVars = cellfun(@numval, problem.equations)';
           cumVars = cumsum(numVars);
           ii = [[1;cumVars(1:end-1)+1], cumVars];
           
           eqn = size(ii,1);
           dx = cell(eqn,1);
           for i = 1:eqn
               dx{i} = result(ii(i,1):ii(i,2));
           end
       end
       
       function solver = setupSolver(solver, A, b, varargin) %#ok 
           % Run setup on a solver for a given system
           
           % Dummy function run before a set of linear problems with
           % different right hand sides
       end
       
       function  solver = cleanupSolver(solver, A, b, varargin) %#ok 
           % Clean up solver after use (if needed)
           
           % Dummy function run after a set of linear problems with
           % different right hand sides 
           
           % For overloading when for example calling multigrid solvers as
           % a preconditioner
       end
   end
end