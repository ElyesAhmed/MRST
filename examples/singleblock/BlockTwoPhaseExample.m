%% Block One-Phase Example
% 
% Example showing permeability two-phase upscaling of a single grid block,
% i.e., the entire grid G is upscaled to a single coarse cell.
% 
% Step through the blocks of the script.
% 

%% Add MRST modules

% We rely on the following MRST modules
mrstModule add incomp upscaling ad-props ad-core

% In addition, we will use the SPE10 grid as an example model
mrstModule add spe10


%% Construct a model
% We extract a small block from the SPE10 model 1.

% Rock
I = 1:5; J = 30:35; K = 1:5; % Make some selection
rock = SPE10_rock(I, J, K); % Get SPE10 rock

% Grid
cellsize = [20, 10, 2].*ft; % Cell size (ft -> m)
celldim  = [numel(I), numel(J), numel(K)];
physdim  = celldim.*cellsize;
G = cartGrid(celldim, physdim); % Create grid structre
G = computeGeometry(G); % Compute volumes, centroids, etc.

% Fluid
fprop = getExampleFluidProps();
fluid = initADIFluidOW(fprop);

% Close figure if open
close(figure(1)); 

%% Plot fluid properties

figure(1); clf; hold on;
sW = linspace(0.1,0.9,100)';
plot(sW, fluid.krW(sW), 'b');
plot(sW, fluid.krO(1-sW), 'r');

%%

figure(1); clf; hold on;
sW = linspace(0.1,0.9,100)';
plot(sW, fluid.pcOW(sW)./barsa, 'b');


%% Plot the grid

% Create a new figure and set it wide
close(figure(1)); fh = figure(1);
op = get(fh, 'OuterPosition');
set(fh, 'OuterPosition', op.*[1 1 2.4 1]);

% Permeability on a log-scale
subplot(1,2,1);
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)) );
view(3); axis tight;
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
title('Permeability'); colorbar;

% Porosity
subplot(1,2,2);
plotCellData(G, rock.poro); view(3); axis tight;
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
title('Porosity'); colorbar;


%% Upscale block using steady-state capillary limit

% Create GridBlock: to simplify the passing of arguments in different
% upscaling methods, we create a GridBlock instance which holds the grid,
% the rock, and other properties of the grid block that is to be upscaled.
block = GridBlock(G, rock);

% Upscale permeability: we first need to upscale the permeability. We here
% use Dirichlet boundary conditions with no-flow on the normal boundaries
KU = upAbsPermPres(block) %#ok<NOPTS>






%% Upscale block using periodic BC

% Create GridBlock with the parameter 'periodic' set to 'true'. The
% GridBlock will then make the grid periodic.
block = GridBlock(G, rock, 'periodic', true);

% Upscale permeability. The function is set up to only return the diagonal
% of the upscaled tensor by default.
KUp = upAbsPermPres(block) %#ok<NOPTS>

% We can get the full tensor by asking for it.
upAbsPermPres(block, 'fulltensor', true)

% Note that the porosity will not change becuase the grid is periodic.


%% Using the Upscaler class

% Instead of calling the upscaling functions directly, we may use the
% subclasses of the Upscaling class. For one phase upscaling, we call the
% class OnePhaseUpscaler.

% Create an instance of the upscaler 
upscaler = OnePhaseUpscaler(G, rock);

% Perform the upscaling. The data structure returned contains both the
% permeability and the porosity.
data = upscaler.upscaleBlock(block) %#ok<NASGU,NOPTS>


%% Permeability averaging

% We may also choose to use an averaging method to find uspcaled values of
% the permeability instead of the pressure solver. Let us for example
% compute the arithmetic average:
upscaler.OnePhaseMethod = 'arithmetic';
data = upscaler.upscaleBlock(block) %#ok<NASGU,NOPTS>

% Another alternative is the combination of harmonic and arithmetic. For
% each different method, observe how the values of the upscaled
% permeability changes.
upscaler.OnePhaseMethod = 'harmonic-arithmetic';
data = upscaler.upscaleBlock(block) %#ok<NOPTS>



