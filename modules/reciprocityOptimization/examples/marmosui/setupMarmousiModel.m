% setupModel: construct Marmousi model (100×100) from the three stored realizations
%   modelChoice: 'first', 'second', 'third', or 'mean' (default: 'mean')
function [G, rock, fluid, W, schedule] = setupMarmousiModel(modelChoice)

    if nargin < 1
        modelChoice = 'mean';   % use average of the three blocks by default
    end

    % --- 1. Grid and geometry definition (100×100) ---
    nx = 100;   % true horizontal cells
    nz = 100;   % vertical cells
    Lx = 9200;  % physical length [m]
    Lz = 3000;  % physical depth [m]
    G = cartGrid([nx, 1, nz], [Lx, 10, Lz]);
    G = computeGeometry(G);

    % --- 2. Read binary Marmousi files (100×300 total) ---
    vfile = 'Marmousi-Velocity-VNXZ300100.dat';
    dfile = 'Marmousi-Density-VNXZ300100.dat';

    % Velocity
    fid = fopen(vfile, 'r');
    if fid == -1, error('Cannot open %s', vfile); end
    fseek(fid, 0, 'eof'); fileBytes = ftell(fid); fseek(fid, 0, 'bof');
    precision = 'float32';
    if fileBytes/(300*100) > 6, precision = 'float64'; end
    v_raw = fread(fid, 300*100, precision);
    fclose(fid);

    % Density
    fid = fopen(dfile, 'r');
    if fid == -1, error('Cannot open %s', dfile); end
    rho_raw = fread(fid, 300*100, precision);
    fclose(fid);

    % Reshape to [nz × nx_total]  => 100 × 300
    v_full   = reshape(v_raw,   [300, 100])';   % 100×300
    rho_full = reshape(rho_raw, [300, 100])';

    % --- 3. Extract the three 100×100 blocks ---
    v1 = v_full(:, 1:100);
    v2 = v_full(:, 101:200);
    v3 = v_full(:, 201:300);

    rho1 = rho_full(:, 1:100);
    rho2 = rho_full(:, 101:200);
    rho3 = rho_full(:, 201:300);

    % --- 4. Choose the desired model ---
    switch lower(modelChoice)
        case 'first'
            v = v1;   rho = rho1;
        case 'second'
            v = v2;   rho = rho2;
        case 'third'
            v = v3;   rho = rho3;
        case 'mean'
            v = (v1 + v2 + v3) / 3;
            rho = (rho1 + rho2 + rho3) / 3;
        otherwise
            error('modelChoice must be ''first'', ''second'', ''third'', or ''mean''.');
    end

    % --- 5. Rock physics transforms (porosity & permeability) ---
    v_ma = 5500;   v_f = 1500;
    phi  = (v_ma - v) ./ (v_ma - v_f);
    phi  = max(0.01, min(phi, 0.4));   % floor at 1% porosity

    grain_d = 100e-6;
    c_perm  = grain_d^2 / 180;
    k = c_perm * (phi.^3) ./ ((1 - phi).^2 + 1e-12);

    % Build MRST rock structure (MRST column‑major ordering)
    rock = makeRock(G, k(:), phi(:));
    rock.poro = phi(:);

    % --- 6. Fluid definition ---
    pRef = 200*barsa;
    fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                               'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                               'n',     [2, 2, 0]);
    c_comp = 1e-5/barsa;
    fluid.bO = @(p) exp((p - pRef)*c_comp);

    % --- 7. Well setup (scaled for 100×100 grid) ---
% --- 7. Realistic two‑well setup (left injector, right producer) ---
% Grid: 100 (x) × 1 × 100 (z), 9200 m × 3000 m → each cell ~92 m × 30 m.
% Injector: left side, completed in deeper water‑leg (cells 80–98, ~2400‑2940 m)
inj_x   = 12;                % left side
inj_z   = 80:98;             % deep completion
inj_cells = sub2ind([nx, 1, nz], repmat(inj_x, 1, numel(inj_z)), ...
                    ones(1, numel(inj_z)), inj_z);

% Producer: right side, completed in shallower oil‑zone (cells 45–70, ~1350‑2100 m)
prod_x  = 88;                % right side
prod_z  = 45:70;             % mid‑depth completion
prod_cells = sub2ind([nx, 1, nz], repmat(prod_x, 1, numel(prod_z)), ...
                     ones(1, numel(prod_z)), prod_z);

W = [];
W = addWell(W, G, rock, inj_cells,  'Type', 'rate', ...
            'Val', 5000*meter^3/day, 'comp_i', [1 0], ...
            'Name', 'I', 'Sign', 1);
W = addWell(W, G, rock, prod_cells, 'Type', 'bhp', ...
            'Val', 150*barsa, 'comp_i', [0 1], ...
            'Name', 'P', 'Sign', -1);

    % --- 8. Simulation schedule ---
% --- 8. Realistic Simulation Schedule (approx. 12.3 years) ---
% Early time: smaller steps to capture initial pressure/saturation changes
ts1 = repmat(30*day, 20, 1);            % 20 x 30 d = 600 days
% Main production phase: moderate steps for flood progression
ts2 = repmat(60*day, 50, 1);            % 50 x 60 d = 3000 days
% Late time: larger steps as system matures
ts3 = repmat(90*day, 10, 1);            % 10 x 90 d = 900 days
% Total time: 4500 days (~12.3 years)

numCnt = 3;   % three control periods, all using the same well setup
[schedule.control(1:numCnt).W] = deal(W);

% All controls identical (constant rate / constant BHP) – you can vary them later
schedule.step.control = rldecode((1:numCnt)', [numel(ts1); numel(ts2); numel(ts3)]);
schedule.step.val     = vertcat(ts1, ts2, ts3);
    %     ts = { [1 1 3 5 5 10 10 10 15 15 15 15 15 15 15]'*day, ...
%            repmat(150/10, 10, 1)*day, ...
%            repmat(150/6, 6, 1)*day, ...
%            repmat(150/6, 6, 1)*day };
%     numCnt = numel(ts);
%     [schedule.control(1:numCnt).W] = deal(W);
%     schedule.step.control = rldecode((1:4)', cellfun(@numel, ts));
%     schedule.step.val     = vertcat(ts{:});
    gravity off;

    disp(['Marmousi model (100x100) loaded using block: ' modelChoice]);
end