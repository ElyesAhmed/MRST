function operators = setMpsaOperators(model)
    
    G = model.G;
    mech = model.mech;
    prop = mech.prop;
    loadstruct = mech.loadstruct;
    
    eta = 0;
    bcetazero = false;
    
    [tbls, mappings] = setupStandardTables(G);

    assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, 'bcetazero', bcetazero, 'adoperators', true); 
    
    operators = assembly.adoperators;
    
end

