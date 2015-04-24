%% Figure 1-5

% These figures are not generated by Matlab

%% Figure 6

plotVEReconstruction_new(6);

%% Figure 7, 9 and 10

plotVEProfilesAndRelperms;

%% Figure 8

plotVEReconstruction_new(8);

%% Figure 11 and 12

% NB: the code for the runStandardModel function is located in the 'paper2' folder.
runStandardModel('data/dissolutionTopSurfaceExample' , ...
                 @plotTopSurfaceExample              , ...
                 'A' , [0 2]                         , ...
                 'residual', true                    , ...
                 'dis_types', {'none', 'rate'}       , ...
                 'fluid_types', {'sharp interface'   , ...
                                 'linear cap.'       , ...
                                 'P-scaled table'    , ...
                                 'P-K-scaled table'});