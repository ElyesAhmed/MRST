function fluid = makeVEFluidsForTest(fluid,fluid_case,varargin)
    opt=struct('res_oil',0,'res_gas',0,'Gt',[],'rock',[]);
    opt = merge_options(opt, varargin{:});
    fluid_names={'simple',...
                 'sharp_interface',...                
                 'cap_linear',...                 
                 'cap_1D_table_P',...
                 'cap_1D_table_kscaled',...
                 'cap_1D_table_SH',...
                 'integrated'};
%     fluid_names={ 'integrated'};       
%             ,...                 
%                 'VE_1D_table_test'                                     
%                 };
%'Uniform_table',...
     if(nargin==0)
         fluid=fluid_names;
         return;
     else
        if(~any(strcmp(fluid_case,fluid_names)))
           disp(['Error wrong fluid name', fluid_case])
           disp('Valid names are')
           for i=1:numel(fluid_names)
              disp(fluid_names{i}); 
           end
           error('Wrong fluid name')
        end
     end


switch fluid_case
    case 'simple'
       fluid.krG=@(sg,varargin) sg;
       fluid.krOG=@(so,varargin) so;
       fluid.pcOG=@(sg, p, varargin) norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*(sg).*opt.Gt.cells.H;       
       fluid.res_gas=0;
       fluid.res_oil=0;
       fluid.invPc3D=@(p) 1-(sign(p+eps)+1)/2;
       fluid.kr3D=@(s) s;
    case 'integrated'
        fluid = addVERelpermIntegratedFluid(fluid, 'res_oil',opt.res_oil,...
                                            'res_gas',opt.res_gas,...
                                            'Gt',opt.Gt,...
                                            'kr_pressure',true,...
                                            'Gt',opt.Gt,...
                                            'int_poro',false,...
                                            'rock',opt.rock);
    
    case 'sharp_interface'    
       fluid = addVERelperm(fluid, opt.Gt, ...
                            'res_oil',opt.res_oil,...
                            'res_gas',opt.res_gas);
    
    %case 'Uniform_table'
    %    pc3D =@(s) 10*barsa./s.^2;
    %    [Sdp_table, dpkrS_table] = makeSdp_table('pc3D',pc3D,'maxdP',H*drho*10*(4));
    %    fluid = addVERelpermTable('Sdp_table',Sdp_table,'dpkrS_table',dpkrS_table,'Gt',Gt);
    case 'cap_linear'
        fluid = addVERelpermCapLinear(fluid,...
                                      'res_gas',opt.res_gas,...
                                      'res_oil',opt.res_oil,...
                                      'beta',2,...
                                      'cap_scale',0.1*max(opt.Gt.cells.H)*10*(fluid.rhoOS-fluid.rhoGS),...
                                      'H',opt.Gt.cells.H,'kr_pressure',true);                                 
        %fluid = addVERelpermCapLinear(fluid,'res_gas',0.1,'beta',4,'cap_scale',0.3*H*10*(fluid.rhoOS-fluid.rhoGS),'H',Gt.cells.H,'kr_pressure',false);
    case 'VE_1D_table_test'
        %
        S_tab=linspace(0,1,10)';
        kr_tab=S_tab;
        h_tab=S_tab*max(opt.Gt.cells.H);
        %h_mean=mean(Gt.cells.H);
        table_co2_1d=struct('SH',S_tab.*opt.Gt.cells.H,'krH',kr_tab.*opt.Gt.cells.H,'h',h_tab,'is_kscaled',false);
        table_co2_1d.invPc3D =@(p) 1-(sign(p+eps)+1)/2;
        table_co2_1d.kr3D=@(s) s;
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTables(fluid,...
                                   'height',opt.Gt.cells.H,...
                                   'table_co2',table_co2_1d,...
                                   'table_water',table_water_1d); 
     case 'cap_1D_table_SH'
         drho=400;
        C=max(opt.Gt.cells.H)*0.4*drho*norm(gravity);
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d=makeVEtables('invPc3D',@(p) max((C./(p+C)).^(1/alpha),opt.res_oil),...
                                  'is_kscaled',false,...
                                  'kr3D',@(s) s.^beta,...
                                   'drho',drho,...
                                   'Gt',opt.Gt,...
                                   'samples',samples);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTables(fluid,...
                                   'res_oil',opt.res_oil,...
                                   'res_gas',opt.res_gas,...
                                   'height',opt.Gt.cells.H,...
                                   'table_co2',table_co2_1d,...
                                   'table_water',table_water_1d);
     case 'cap_1D_table_P'
        drho=400;
        C=max(opt.Gt.cells.H)*0.4*drho*norm(gravity);
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('invPc3D', @(p) max((C./(p+C)).^(1/alpha),opt.res_oil),...
            'is_kscaled', false,'kr3D', @(s) s.^beta,...
            'drho', drho,...
            'Gt', opt.Gt, 'samples', samples);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,...
                                           'res_oil', opt.res_oil,...
                                           'res_gas', opt.res_gas,...
                                           'height',opt.Gt.cells.H,...
                                           'table_co2',table_co2_1d,...
                                           'table_water',table_water_1d,...
                                           'kr_pressure',true); 
      case 'cap_1D_table_kscaled'        
        %surface_tension=100;  
        kscale=sqrt(0.1/(100*milli*darcy))*fluid.surface_tension;
        drho=400;
        %C=max(opt.Gt.cells.H)*0.4*drho*norm(gravity)/kscale;
        C=1;
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('invPc3D', @(p) max((C./(p+C)).^(1/alpha),opt.res_oil),...
                                    'is_kscaled', true,....
                                    'kr3D', @(s) s.^beta,...
                                    'drho', drho,...
                                    'Gt', opt.Gt,...
                                    'samples', samples,'kscale',kscale);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,...
                                           'res_oil', opt.res_oil,...
                                           'res_gas', opt.res_gas,...
                                           'height', opt.Gt.cells.H,...
                                           'table_co2',table_co2_1d,...
                                           'table_water',table_water_1d,...
                                           'rock',opt.rock,...
                                           'kr_pressure',true); 
    otherwise
       error('No such fluid case')
end
fluid.name=fluid_case;