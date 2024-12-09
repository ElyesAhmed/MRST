clear; clc;
mrstModule add biochemistry compositional ad-blackoil ad-core ad-props mrst-gui
gravity reset on
useHandler=false;

%====Geometrie, maillage et proprietes roche
[nx,ny,nz] = deal(20,1,20);
[Lx,H]=deal(200,100);
dims = [nx, ny,nz];
pdims = [Lx, 1, H];
dz=H/nz;
G = cartGrid(dims, pdims); 
depth_res=1100; %profondeur haut reservoir
G.nodes.coords(:,3) = G.nodes.coords(:,3) + depth_res;

G = computeGeometry(G);
rock = makeRock(G, 30*milli*darcy, 0.2);

[i1,i2,i3,i4]=deal(floor(0.05*nz),floor(0.4*nz),floor(0.5*nz),floor(0.7*nz));
for i=1:i1
    I = (1:nx).'; 
    J = ones(nx,1); %repmat(1, [nx, 1]);
    K = repmat(i, [nx, 1]);
    ind = sub2ind(G.cartDims, I, J, K);
    rock.perm(ind) = 1e-19;
    rock.poro(ind) = 0.001;
end
for i=i2+1:i3
    I = (1:nx).'; 
    J = ones(nx,1);%repmat(1, [nx, 1]);
    K = repmat(i, [nx, 1]);
    ind = sub2ind(G.cartDims, I, J, K);
    rock.perm(ind) = 1e-15;
    rock.poro(ind) = 0.05;
end
 figure;
 plotCellData(G,rock.poro,'EdgeAlpha',0.2)
 colorbar
 view(3)


%%=========initialisation proprietes du fluide
%===Compositional fluid model (initialization with the CoolProp library)=
compFluid =TableCompositionalMixture({'Water','Hydrogen','CarbonDioxide',...
    'Nitrogen','Methane'}, {'H2O','H2','CO2','N2','CH4'});
disp (compFluid)

%Liquid-Gas (H2)
[rhol,rhog]=deal(1000* kilogram/meter^3,8.1688* kilogram/meter^3); %density kilogram/meter^3;
[viscol,viscog]=deal(1.0*centi*poise,0.0094234*centi*poise);%viscosity
[cfl,cfg]=deal(4.521e-5/barsa,8.1533e-3/barsa); %compressibility
 
[srl,src]=deal(0.0,0.0);
Phydro0=rhol*norm(gravity).*G.cells.centroids(:,3);
[Pmaxz,Pref1,Pminz,Pe]=deal(95*barsa,114*barsa,120*barsa,0.1*barsa); %pressions

%initialisation fluides incompressibles, Brooks-Corey relperm krw=(Sw)^nw
fluid=initSimpleADIFluid('phases', 'WG', 'mu',[viscol,viscog],...
    'rho',[rhol,rhog],'pRef',Pref1,'c',[0.0,cfg],'n',[2,2],'smin',[srl,src]);
% Pression capillaire
pcWG = @(sw) Pe * sw.^(-1/2);
fluid.pcWG = @(sg) pcWG(max((1-sg-srl)./(1-srl), 1e-5)); 


T = 100*day;
pv=sum(poreVolume(G,rock))/T;
rate = 100*pv;
niter = 100;%300

%%==============Conditions aux limites et puits============
bc = [];
%puit d'injection
W = [];
W = verticalWell(W, G, rock,1,1,nz, 'comp_i', [0, 1],'Radius',0.5,...
    'name', 'Injector', 'type', 'rate','Val',rate, 'sign', 1);
W(1).components = [0.0 0.95 0.05 0.0 0.0];

%% =======Model Setup: Compositional Model with Bacterial Growth=======
arg = {G, rock, fluid, compFluid,'water', true, 'oil', false, ...
    'gas', true, 'bacteriamodel', true,'bDiffusionEffect',false,...
    'liquidPhase', 'W','vaporPhase', 'G'}; % water=liquid, gas=vapor
model = BiochemistryModel(arg{:});
model.outputFluxes = false;
eosname='pr'; %'sw';
model.EOSModel = SoreideWhitsonEquationOfStateModel(G, compFluid,eosname);
model.EOSModel.msalt=0;
            
%===Conditions initiales=====================
T0=317.5;
s0= [0.8 0.2]; %initial saturations  Sw=1
z0 = [0.8,0.0,0.006,0.018,0.176]; %initial composition: H2O,H2,CO2,N2,CH4.

%===Bacteria model===========================
nbact0=10^6;
state0 = initCompositionalStateBacteria(model,Phydro0,T0,s0,z0,nbact0);

%===Resolution pression/transport========================
deltaT = T/niter;
schedule = simpleSchedule(repmat(deltaT,1,niter),'bc', bc,'src', [],'W',W);
nls = NonLinearSolver('useRelaxation', true);


if useHandler
    dir='/home/sdelage2/PROJETS/gdr_h2/MRST2024_microbes/MRST/output';
    diroutput='WaterGaz_H2OH2';%'OilGaz_H2OH2'
    handler = ResultHandler('writeToDisk', true,'dataDirectory',dir,...
        'dataFolder', diroutput);
    [~,states,report] = simulateScheduleAD(state0, model, schedule,...
        'nonlinearsolver', nls,'outputHandler', handler);
else
    [~,states,report] = simulateScheduleAD(state0, model, schedule,'nonlinearsolver', nls);
end  

%===Plottings============================================
if useHandler
    %Extraction data states de handler
    handler1 = ResultHandler('dataDirectory',dir,'dataFolder','WaterGaz_H2OH2');
    m = handler1.numelData();
    states = cell(m, 1);
    for i = 1:m
        states{i} = handler1{i};
    end
end

time=0;
figure;
for i= 1:niter
    x = G.cells.centroids(:,1);
    z = G.cells.centroids(:,3);
    X = reshape(x, [nx,nz]);
    Z = reshape(z, [nx,nz]);
    xH2 = reshape(states{i}.x(:,2), [nx,nz]); %in liquid phase
    yH2 = reshape(states{i}.y(:,2), [nx,nz]); %in gaz phase 
    zH2 = reshape(states{i}.components(:,2), [nx,nz]);%overall
    zH2O = reshape(states{i}.components(:,1), [nx,nz]);
    Sw = reshape(states{i}.s(:,1), [nx,nz]);
    Sg = reshape(states{i}.s(:,2), [nx,nz]);
    Pres= reshape(states{i}.pressure, [nx,nz]);
    nbacteria=reshape(states{i}.nbact, [nx,nz]);

    subplot(2,3,1);   
    contourf(X,Z,Sw,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Water saturation','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,3,4);   
    contourf(X,Z,Sg,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Gaz saturation','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,3,2);   
    contourf(X,Z,zH2,60,'EdgeColor','auto');
    clim([0 0.8])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('Overall H2 molar fraction','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar 

    subplot(2,3,3); 
    contourf(X,Z,nbacteria,60,'EdgeColor','auto');
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('microbial density','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar

    subplot(2,3,5); 
    contourf(X,Z,xH2,60,'EdgeColor','auto');
    clim([0 7.e-5])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('H2 molar fraction in liquid','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar
   
    subplot(2,3,6); 
    contourf(X,Z,yH2,60,'EdgeColor','auto');
    clim([0 1])
    axis equal
    axis ([0 Lx  depth_res depth_res+H])
    ylabel('H2 molar fraction in gaz','FontSize',15)
    set(gca,'Ydir','reverse')
    colormap('jet')
    colorbar
  
    time = time + deltaT;
    title(sprintf('injection duration = %.2f days',convertTo(time,day)))
    pause(0.001)
end