function obj=h2oProps_simple()

obj.density =@(p,T) 1000+0*p;
obj.enthalpy =@(p,T)  4180*(T - 293.15);
obj.viscosity = @(p,T)  1e-03+0*p;
obj.name='h2o simple';
end