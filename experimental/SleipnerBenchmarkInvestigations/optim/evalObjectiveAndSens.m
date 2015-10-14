function [val, der, wellSols, states] = evalObjectiveAndSens(u, obj, state0, model, schedule, scaling)
% Objective (and gradient) evaluation function based on input control vector u
% u(1:n) : dz
% u(n+1) : rhomult
% u(n+2) : permmult
% u(n+3)  : poromult

% where minimizing -> set objective to minus objective:

minu = min(u);
maxu = max(u);
if or(minu < -eps , maxu > 1+eps)
    warning('Controls are expected to lie in [0 1]\n')
end

boxLims = scaling.boxLims;
if isfield(scaling, 'obj')
    objScaling = scaling.obj;
else
    objScaling = 1;
end

% update model
% dz, rhofac, permfac, porofac

n = model.G.cells.num;
% 1.

% multipliers
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
us = u.*(umax-umin)+umin;
for i=1:numel(schedule.control)
    schedule.control(i).dz       = us(1:n);
    schedule.control(i).rhofac   = us(n+1);
    schedule.control(i).permfac  = us(n+2);
    schedule.control(i).porofac  = us(n+3);
end

% run simulation:
[wellSols, states, sim_report] = simulateScheduleAD(state0, model, schedule);

states = addHeightData(states, model.G, model.fluid);
% compute objective:
vals = obj(wellSols, states, schedule);
val  = - sum(cell2mat(vals))/objScaling;

% run adjoint:
if nargout > 1
    objh = @(tstep)obj(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
    g    = computeGradientAdjointAD(state0, states, model, schedule, objh, 'ControlVariables', {'scell','mult'});
    g    = cell2mat(g);
    g    = sum(g,2);
    % scale gradient:
    dBox   = boxLims(:,2) - boxLims(:,1);
    der  = - (dBox/objScaling).*g;
end
end

% function grd = scaleGradient(grd, schedule, boxLims, objScaling)
% dBox   = boxLims(:,2) - boxLims(:,1);
% for k = 1:numel(schedule.control)
%     grd{k} = (dBox/objScaling).*grd{k};
% end
% end
    
