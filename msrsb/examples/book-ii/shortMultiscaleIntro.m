%% Short introduction to multiscale finite-volume methods
% This example is an abbreviated version of the introMultiscale.m script
% that has been edited so that it only contains the code necessary to
% produce the plots presented in Section 2.5 of the multiscale tutorial
% chapter.

mrstModule add coarsegrid spe10 msrsb incomp upscaling

%% Setup and solve fine-scale model
% As an example of a strongly heterogeneous problem, we consider a 20x20
% subsample from Model 2 of the SPE10 upscaling benchmark with a pressure
% drop from the west to the east boundar.
G      = computeGeometry(cartGrid([20 20], [20 20].*[20 10]*ft));
rock   = getSPE10rock(1:20, 1:20, 15);
rock.perm = rock.perm(:, 1);
hT     = computeTrans(G, rock);
bc     = pside([], G, 'West', 100*barsa);
bc     = pside(bc, G, 'East',  50*barsa);
fluid  = initSingleFluid('rho', 1, 'mu', 1);
state0 = initResSol(G, 0);
state  = incompTPFA(state0, G, hT, fluid, 'MatrixOutput', true, 'bc', bc);
n      = G.cells.num;
A      = full(abs(state.A));

%% Multiscale solution
cdims = [5, 4];
CG    = coarsenGeometry(generateCoarseGrid(G, partitionUI(G, cdims)));
CG    = storeInteractionRegion(CG);
basis = getMultiscaleBasis(CG, A);
Am    = basis.R*state.A*basis.B;
bm    = basis.R*state.rhs;
p_ms  = Am\bm;
m     = CG.cells.num;

%% Upscaled solution
crock.perm = upscalePerm(G, CG, rock);
chT    = computeTrans(CG, crock);
cbc    = coarsenBC(CG, bc);
cstate = initResSol(CG, 0);
cstate = incompTPFA(cstate, CG, chT, fluid, 'MatrixOutput', true, 'bc', cbc);

%% Compare the solutions
figure('Position',[300 340 1000 420]);
ca = [min(state.pressure), max(state.pressure)]/barsa;
colormap default
subplot(2, 2, 1);
plotCellData(CG, cstate.pressure/barsa,'EdgeColor','none'); 
plotFaces(CG,1:CG.faces.num)
axis equal tight; caxis(ca); title('Upscaled')

subplot(2, 2, 2);
plotCellData(G, state.pressure/barsa);
axis equal tight; caxis(ca); title('Fine-scale')

subplot(2, 2, 3);
plotCellData(CG, p_ms/barsa,'EdgeColor','none'); 
plotFaces(CG,1:CG.faces.num)
axis equal tight; caxis(ca); title('Coarse MS');

subplot(2, 2, 4);
plotCellData(G, basis.B*p_ms/barsa);
axis equal tight; caxis(ca); title('Fine MS');
set(colorbar,'Position',[.9 .11 .025 .34]);
colormap(parula(10))

%% Computed versus reconstructed velocity field
% Compute velocity field directly from multiscale pressure solution on fine
% scale and compare with a reconstructed version. As a measure, we use the
% discrete divergence.
mrstModule add ad-core
rec  = incompMultiscale(state0, CG, hT, fluid, basis, 'bc', bc, 'reconstruct', true);
prol = incompMultiscale(state0, CG, hT, fluid, basis, 'bc', bc, 'reconstruct', false);

op   = setupOperatorsTPFA(G, rock);
Div = @(flux) op.Div(flux(op.internalConn));
[df, drec, dprol] = deal( Div(state.flux), Div(rec.flux), Div(prol.flux) );

figure('Position',[80 340 1400 420]);
og = generateCoarseGrid(G,ones(G.cells.num,1));
for i=1:3
    switch i
        case 1
            div = Div(state.flux); tittel='Fine-scale flux';
        case 2
            div = Div(rec.flux); tittel='Reconstructed flux';
        case 3
            div = Div(prol.flux); tittel='Flux from prolongated pressure';
    end
    subplot(1,3,i)
    plotCellData(G, div, 'EdgeColor','none'); axis equal tight
    caxis([-.46e-7 .46e-7]); title(tittel);
    plotFaces(og,1:og.faces.num,'LineWidth',1);
end
cmap = interp1([0; 0.5; 1], [1, 0, 0; 1, 1, 1; 0, 0, 1], 0:0.01:1);
colormap(cmap);
set(colorbar,'Position',[.925 .33 .015 .37]);

%% Plot fine-scale pressures side by side - together with streamlines
% We have excellent agreement between the fine-scale pressures and the
% velocity field when visualized as streamlines.
mrstModule add streamlines
ijk   = gridLogicalIndices(G);
cells = find(ijk{1} == 1);
p     = repmat([0.25, 0.25; 0.75, 0.25; ...
                0.25, 0.75; 0.75, 0.75; 0.5, 0.5], numel(cells), 1);
trace = @(state) pollock(G, state, [rldecode(cells,5), p], 'pvol', rock.poro);


figure('Position',[300 340 1000 420]);
subplot(1, 2, 1);
plotFaces(og,1:og.faces.num,'LineWidth',1); axis equal tight
h = streamline(trace(state));
set(h, 'Color', 'k','linewidth',.5); title('Fine-scale');

subplot(1, 2, 2);
plotFaces(og,1:og.faces.num,'LineWidth',1); axis equal tight
h = streamline(trace(rec));
set(h, 'Color', 'k','linewidth',.5); title('Multiscale');

%% Compare MsFV and MsRSB
figure('Position',[300 340 1000 420]);
subplot(1,2,1)
err = @(ms) (ms.pressure - state.pressure)./state.pressure;
plotCellData(G, err(rec),'EdgeColor','none'); axis equal tight
plotFaces(og,1:og.faces.num,'LineWidth',1);
caxis([-.151 .151]), title('MsRSB');

subplot(1,2,2)
wb      = partitionUIdual(CG,[5 4]);
CG.dual = makeExplicitDual(CG, wb);
basis2  = getMultiscaleBasis(CG, A, 'type', 'msfv');
msfv    = incompMultiscale(state0, CG, hT, fluid, basis2, 'bc', bc);
plotCellData(G, err(msfv),'EdgeColor','none'); axis equal tight
plotFaces(og,1:og.faces.num,'LineWidth',1);
caxis([-.151 .151]), title('MsFV');
colormap(cmap);
set(colorbar,'Position',[.925 .32 .015 .39]);