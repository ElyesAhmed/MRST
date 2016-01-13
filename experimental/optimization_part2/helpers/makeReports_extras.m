function reports = makeReports_extras(Gt, states, rock, fluid, schedule, ...
                               residual, traps, dh, wellSols)
% residual on form: [sw, sr]
% states on form:   {initState, states{:}} to include initial state

   assert( numel(states) == numel(schedule.step.val)+1 , ...
       'Ensure the initial state has been included in the varargin ''states''.')
   
   tot_inj = 0;
   for i = 1:numel(states)

      [h, h_max] = computePlumeHeight(Gt, states{i}, residual(1), residual(2));

      if i == 1
         reports(i).t         = 0;
         reports(i).W         = [];
      else
         reports(i).t = sum(schedule.step.val(1:i-1));
         reports(i).W = schedule.control(schedule.step.control(i-1)).W;
         
         %assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         %tot_inj = tot_inj + sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS;
         
         % since bhp-controlled wells may have been used, we cannot rely on
         % using W.val for the rate. Thus, we use the SURFACE RATE found in
         % the wellSols structure.
         % NB: for rate-controlled wells, wellSols.val == wellSols.qGs, not wellSols.qGr !!!
         tot_inj = tot_inj + sum([wellSols{i-1}.qGs]) * schedule.step.val(i-1) * fluid.rhoGS;
         
      end
      
      reports(i).sol       = states{i};
      reports(i).sol.h     = h;
      reports(i).sol.h_max = h_max;
      
      reports(i).masses    = massTrappingDistributionVEADI(Gt          , ...
                                                        reports(i).sol , ...
                                                        rock           , ...
                                                        fluid          , ...
                                                        traps          , ...
                                                        dh);
      leaked = tot_inj - sum(reports(i).masses);
      reports(i).masses = [reports(i).masses, leaked];

   end
   
   % hack @@ put well cell indexes in reports(1) using reports(2)
      reports(1).W = reports(2).W;
   
end

