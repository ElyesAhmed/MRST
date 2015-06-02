G = computeGeometry(cartGrid([8 8 8]));
coarse = [1 1 1];
cG = computeGeometry(cartGrid(coarse,G.cartDims));
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);
fl = {'East', 'North', 'Top'};
fr = {'West', 'South', 'Bottom'};
clc
disp('                      Arithmetic        Harmonic      Harm-arith');
for n=1:3
   disp('-----------------------------------------------------------------');
   switch n
      case 1
         K = repmat([50 150 400 300 10 250 100 350], ...
            prod(G.cartDims(1:2)),[]);
         rock.perm = K(:)*[1 1 1]*milli*darcy;
         disp('Layered:');
      case 2
         disp('Tarbert');
         rock = SPE10_rock(9:16,9:16,3:10);
         rock.perm = rock.perm.*milli*darcy;
      case 3
         disp('Upper Ness');
         rock = SPE10_rock(9:16,9:16,43:50);
         rock.perm = rock.perm.*milli*darcy;
   end

   subplot(1,3,n);
   plotCellData(G,log10(rock.perm(:,1)/(milli*darcy))); view(3); axis off; zoom(1.3);
   
   %% Upscale values
   q   = ones(G.cells.num,1);
   vol = G.cells.volumes;
   crock = cell(3,1);
   for i=1:size(rock.perm,2)
      % Arithmetic
      crock{1}.perm(:,i) = accumarray(q,vol.*rock.perm(:,i)) ./ ...
         accumarray(q,vol);

      % Harmonic
      crock{2}.perm(:,i) = accumarray(q,vol) ./ ...
         accumarray(q,vol./rock.perm(:,i));

      % Harmonic-arithmetic
      dims = G.cartDims; dims(i)=coarse(i);
      qq = partitionUI(G, dims);
      K = accumarray(qq,vol)./accumarray(qq,vol./rock.perm(:,i));
      crock{3}.perm(:,i) = accumarray(q,K(qq).*vol)./accumarray(q,vol);
   end

   %% Compare upscaling
   for j=1:numel(fl)
      % Fine-scale problem
      bc    = pside([], G, fl{j}, 0);
      faces = bc.face;
      bc    = pside(bc, G, fr{j}, barsa());
      hT    = computeTrans(G, rock);
      xr    = incompTPFA(initResSol(G,0), G, hT, fluid, 'bc', bc);
      flux  = sum(xr.flux(faces));

      fprintf('  %5s -> %6s: ', fl{j}, fr{j});
      for i=1:3
         cbc    = pside([], cG, fl{j}, 0);
         cfaces = cbc.face;
         cbc    = pside(cbc, cG, fr{j}, barsa);
         chT    = computeTrans(cG, crock{i});
         x      = incompTPFA(initResSol(cG,0), cG, chT, fluid, 'bc', cbc);
         fprintf('\t%f', x.flux(cfaces) / flux);
      end
      fprintf('\n');
   end
end
disp('-----------------------------------------------------------------');
cax = caxis;
for i=1:3, subplot(1,3,i), caxis(cax); end
h = colorbar;
set(h,'Position',[0.94 0.1100 0.02 0.8150], ...
   'XTick', .5, 'XTickLabel','[mD]', ...
   'YTickLabel',num2str(10.^(-2:3)'));