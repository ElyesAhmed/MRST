clc; clear; close all;

addpath('../')
addpath('/home/strene/Documents/master/coop/pebiGridding/voronoi3D')

n = 4;
gridLim = [1,1,1];

G = cartGrid([n,n,n],gridLim);
% G = tetrahedronCube([n,n,n], gridLim, 1);
% G = voronoiCube(250  ,gridLim);


% G = computeGeometry(G);
% 
% remCells = (1:G.cells.num)';
% remCells = remCells(G.cells.centroids(:,3) > 1 - ...
%                     (G.cells.centroids(:,1) + G.cells.centroids(:,2))   );
% G = removeCells(G,remCells);
% 
% R1 = @(theta) [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];   
% R2 = @(phi)   [cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
% R3 = @(psi)   [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];  
% 
% X = G.nodes.coords;
% 
% % X(:,2) = X(:,2) + 1/10*X(:,3).^2;
% 
% theta = pi/6; phi = pi/3; psi = pi/18;
% X = X*R1(theta)*R2(phi)*R3(psi);
% 
% G.nodes.coords = X;



%--------------------------------------------------------------------------
%   -\delta u = 0,
%           u = 1/(2\pi||x-C||)
%--------------------------------------------------------------------------
f = @(X) zeros(size(X,1),1);
C = -[.2,.2,.2];
gD = @(X) -1./(2*pi*sqrt(sum((X-repmat(C,size(X,1),1)).^2,2)));
k = 1;

% %--------------------------------------------------------------------------
% %   -\delta u = 1,
% %           u = -(x^2 + y^2 + z^2)/6
% %--------------------------------------------------------------------------
% f = @(X) ones(size(X,1),1);
% gD = @(X) -(X(:,1).^2 + X(:,2).^2 + X(:,3).^2)/6;
% k = 1;


% %--------------------------------------------------------------------------
% %   -\delta u = \sin(x)\cos(y)z(1+alpha^2\pi)(1+\alpha^2\pi^2)  ,
% %           u = \sin(x)\cos(y)z(1+alpha^2\pi)
% %--------------------------------------------------------------------------
% alpha = 2;
% f = @(X) sin(X(:,1)).*cos(alpha*pi*X(:,2)).*X(:,3)*(1+alpha^2*pi^2);
% gD = @(X) sin(X(:,1)).*cos(alpha*pi*X(:,2)).*X(:,3);

% beta = 4*1.8395265*10e-5;
% 
% alpha = gridLim(1)/n*(1/20*beta +1/5)*3*ones(G.cells.num,1);

G = computeVEMGeometry(G,f,k);


alpha = G.cells.diameters*sqrt(3)*(1/20*0.4+1/5);


boundaryFaces = (1:G.faces.num)';
                           
boundaryFaces = boundaryFaces( G.faces.neighbors(:,1) == 0 | ...
                               G.faces.neighbors(:,2) == 0 );

bc = struct('bcFunc', {{gD}}, 'bcFaces', {{boundaryFaces}}, 'bcType', {{'dir'}});

sol = VEM3D(G,f,bc,k);
U = [sol.nodeValues; sol.edgeValues; sol.faceMoments; sol.cellMoments];


Kc = G.cells.centroids;
cells = 1:G.cells.num;
r = .7; c = [1,0,0];
cells = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
faces = G.cells.faces(faceNum);
% faces = 1:G.faces.num;

if k == 2
figure();
plotFaces(G,faces,sol.faceMoments(faces));

colorbar;
view(3);
axis equal;

IF = polygonInt3D(G,1:G.faces.num,gD, 7);
IC = polyhedronInt(G,1:G.cells.num,gD, 7);

u = [gD([G.nodes.coords; G.edges.centroids]); IF./G.faces.areas; IC./G.cells.volumes];

err = abs((U - u));

elseif k == 1
    u = gD(G.nodes.coords);
    err = abs(U-u);
end

h = sum(G.cells.diameters)/G.cells.num;

fprintf('Error: %d\n', norm(err, 2));
figure()
plot(err);  
% % plot(nodeValues)
% 
% %   Implement: efficient rule for faceIntegrals.
% %              change bc to give avg values.




% vols1 = zeros(G.cells.num,1);
% vols2 = zeros(G.cells.num,1);
% 
% 
% V  = [-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
%        0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0); ...
%        1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0); ...
%        0.0,  0.0          ,  3.0/sqrt(6.0)];   
% Vdiff = V(1:end-1,:) - V(2:end,:);
%    
% vol = sqrt(sqrt(8.0)/3.0);
% 
% for i = 1:G.cells.num
%     nodeNum = G.cells.nodePos(i):G.cells.nodePos(i+1)-1;
%     nodes = G.cells.nodes(nodeNum);
%     X = G.nodes.coords(nodes,:);
%     tri = delaunay(X);
%     nTri = size(tri,1);
%     for j = 1:nTri
%         Xdiff = X(tri(j,1:3),:) - X(tri(j,2:4),:);
%         A1 = Vdiff\Xdiff;
%         b = X(tri(j,1),:) - V(1,:)*A1; 
%         Xhat = V*A1 + repmat(b,size(V,1),1);
%         X(tri(j,:),:) - Xhat;
%         
%         A2 = X(tri(j,2:4),:) - repmat(X(tri(j,1),:),3,1);
%         
%         vols1(i) = vols1(i) + abs(det(A1))*vol;
%         vols2(i) = vols2(i) + abs(det(A2))/6;
%         
%     end
% end
% G.cells.volumes - vols1;
% G.cells.volumes - vols2;
% 
% % G.cells.volumes = vols2;
