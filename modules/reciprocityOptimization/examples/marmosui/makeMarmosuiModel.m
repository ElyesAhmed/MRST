%% Load Marmousi from binary .dat files and convert to rock properties (MRST)
% Assumes the files contain a single array (300x100) in row‑major order.


% Grid definition
nx = 300;    % number of cells in x
nz = 100;    % number of cells in z
Lx = 9200;   % physical length [m]
Lz = 3000;   % physical depth [m]

% --- Read binary files ---
vfile = 'Marmousi-Velocity-VNXZ300100.dat';
dfile = 'Marmousi-Density-VNXZ300100.dat';

% Velocity
fid = fopen(vfile, 'r');
fseek(fid, 0, 'eof'); fileBytes = ftell(fid); fseek(fid, 0, 'bof');
bytesPerElement = fileBytes / (nx*nz);
precision = 'float32';  % default to float32; adjust if 240 kB
if bytesPerElement > 6, precision = 'float64'; end
v = fread(fid, nx*nz, precision);
fclose(fid);

% Density (same logic)
fid = fopen(dfile, 'r');
rho = fread(fid, nx*nz, precision);
fclose(fid);

% Reshape: row‑major → [nz, nx] matrix
v   = reshape(v, [nx, nz])';   % now v(nz, nx)
rho = reshape(rho, [nx, nz])';

% Quick check
figure; imagesc(v); colorbar; title('Velocity (raw)');

% --- MRST grid and rock physics ---
G = cartGrid([nx, nz], [Lx, Lz]);
G = computeGeometry(G);

v_ma = 5500; v_f = 1500;
phi  = (v_ma - v) ./ (v_ma - v_f);
phi  = max(0, min(phi, 0.4));

grain_d = 100e-6;
c = grain_d^2 / 180;
k = c * (phi.^3) ./ ((1 - phi).^2 + 1e-12);

rock = makeRock(G, k(:), phi(:));
rock.poro = phi(:);

% Visualise
figure;
plotCellData(G, phi(:)); title('Porosity'); colorbar; axis equal tight;
figure;
plotCellData(G, log10(k(:)./9.869233e-16)); title('log_{10} Permeability (mD)'); colorbar; axis equal tight;
figure;
subplot(1,3,3); plotCellData(G, rho(:)); title('Density'); colorbar; axis equal tight;
colormap(jet);

disp('Done. rock structure ready for simulation.');