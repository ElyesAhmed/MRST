classdef SparseTensor
   
   properties
      % Cell array of components.  Each component should be a struct with 
      % the fields: 
      % - coefs (vector of coefs)
      % - ixs (matrix where each row constitute a multiindex)
      % - indexnames (names of indices, one for each row of 'ixs')
      components 
   end
   
   methods
      
      function self = SparseTensor(varargin)
         switch nargin
           case 1
             % input is directly in the form of a component
             assert(isstruct(varargin{1}));
             self.components = {varargin{1}};
           case 2
             % input is in form of a matrix (or vector) and a list of (one or
             % two) index names
             coefs = varargin{1};
             indexnames = varargin{2}';
             self.components = SparseTensor.make_matrix_tensor(coefs, indexnames);
           case 3
             % input is in the form of coefs, ixs and indexnames.  Coefs may
             % be empty, for indicator tensors
             comp.indexnames = varargin{3};
             comp.ixs = varargin{2};
             comp.coefs = varargin{1};
             comp.coefs = comp.coefs(:); % ensure column vector
             if isempty(comp.coefs)
                comp.coefs = ones(size(comp.ixs, 1), 1);
             end
             self = SparseTensor(comp);
           otherwise
             error('Unsupported arguments to constructor.');
         end
      end
      
      function self = product(self, other, only_semiproduct)
         if nargin < 3
            only_semiproduct = false;
         end

         % if second tensor has existing, contracting indices whose name
         % overlaps with those used in the first tensor, rename them
         [~, cix1] = self.indexNames();
         [~, cix2] = other.indexNames();
         common_cixs = intersect(cix1, cix2);
         for cix = 1:numel(common_cixs)
            new_cix = self.next_unused_contr_ix_name(other);
            other = other.changeIndexName(common_cixs{cix}, new_cix);
         end
         
         % find common index names, and give them new values (since they should be
         % considered to be summed)
         common_names = intersect(self.indexNames(), other.indexNames());
         for nix = 1:numel(common_names)
            % choose an index name that does not overlap with any already
            % used in either of the involved tensors
            flagged_name = self.next_unused_contr_ix_name(other);
            if only_semiproduct
               self = self.duplicateIndex(common_names{nix}, flagged_name);
            else
               self = self.changeIndexName(common_names{nix}, flagged_name);
            end
            other = other.changeIndexName(common_names{nix}, flagged_name);
         end
         self.components = [self.components, other.components];
      end
      
      function self = contract(self, ixname1, ixname2, only_semicontract)
      % can be called with one index (contraction) or two indices
      % (contraction or semicontraction)
         if nargin < 3
            % Simple contraction in one index.  
            comp_ix = self.component_with_ix(ixname1);
            self.components{comp_ix} = ...
                SparseTensor.single_component_contraction(self.components{comp_ix}, ...
                                                         ixname1);
            return
         
         elseif nargin < 4
            % unless requested otherwise, do a full contraction            
            only_semicontract = false;
         end
         
         % contraction (or semicontraction) involving two indices
         newname = self.next_unused_contr_ix_name();
         
         if only_semicontract
            self = self.duplicateIndex(ixname1, newname);
         else
            self = self.changeIndexName(ixname1, newname);
         end
         self = self.changeIndexName(ixname2, newname);
         
         % if the contraction only involves a single component, carry it out
         comp_ix = self.component_with_ix(newname, false);         
         comp = self.components{comp_ix};
         comp = SparseTensor.single_component_contraction(comp);
         self.components{comp_ix} = comp;
      end
      
      function self = plus(self, other)
         self = SparseTensor.apply_binary_operator(self, other, @plus);
      end
      
      function self = minus(self, other)
         self = SparseTensor.apply_binary_operator(self, other, @minus);
      end
      
      function self = rdivide(self, other)
         self = SparseTensor.apply_binary_operator(self, other, @rdivide);
      end
      
      function self = times(self, other)
         self = SparseTensor.apply_binary_operator(self, other, @times);
      end
      
      function self = mpower(self, other)
         self = self.product(other, true);
      end
      
      function t = mtimes(self, other)
         t = self.product(other, false);
      end
      
      function self = sortIndices(self, ixset_order)
      
         self = self.expandall();

         assert(SparseTensor.is_permutation(ixset_order, self.indexNames()));

         perm = SparseTensor.get_permutation(self.indexNames(), ixset_order);

         [self.components{1}.ixs, I] = sortrows(self.components{1}.ixs, perm);

         self.components{1}.coefs = self.components{1}.coefs(I);
         
      end
      
      function self = duplicateIndex(self, name, newname)

         comp_ix = self.component_with_ix(name, false);
         comp = self.components{comp_ix};
         
         ind = strcmp(name, comp.indexnames);
         assert(sum(ind)==1);

         % Add new index values (new column to 'ixs')
         comp.ixs = [comp.ixs, comp.ixs(:,ind)];
         comp.indexnames = [comp.indexnames, newname];
         
         self.components{comp_ix} = comp;
      end
            
      function self = changeIndexName(self, oldnames, newnames)
         if ~iscell(oldnames)
            oldnames = {oldnames}; newnames = {newnames};
         end

         for i = 1:numel(oldnames)
            cur_oldname = oldnames{i};
            cur_newname = newnames{i};

            for c = 1:numel(self.components)
               found = strcmp(cur_oldname, self.components{c}.indexnames);
               if any(found)
                  self.components{c}.indexnames{found} = cur_newname;
               end
            end
         end
      end

      function self = toInd(self)
         self = self.expandall();
         self.components{1}.coefs = self.components{1}.coefs * 0 + 1;
      end
      
      function [ixnames, cixnames] = indexNames(self)
      % return tensor index names.  Contracting indices are returned as
      % second value.
         ixnames = {};
         cixnames = {};
         for c = self.components
            c = c{:};  %#ok
            for name = c.indexnames
               name = name{:}; %#ok
               if SparseTensor.is_contracting_ix(name)
                  cixnames = [cixnames, name]; %#ok
               else
                  % regular index name
                  ixnames = [ixnames, name]; %#ok
               end
            end
         end
      end
      
      function M = asMatrix(self, ixnames, force_sparse, shape)
         SPARSE_THRESHOLD = 200;
         
         if ~exist('force_sparse', 'var')
            force_sparse = false;
         end

         self = self.expandall();

         if ~exist('shape', 'var')
            shape = max(self.components{1}.ixs);
         end
         
         % if the tensor is an intrinsic scalar, print its value
         if numel(self.components{1}.indexnames) == 0
            M = self.components{1}.coefs;
            return
         end
         
         % put input variable 'ixnames' on "standard" form
         if ischar(ixnames)
            ixnames = {{ixnames}};
         else
            assert(iscell(ixnames) && numel(ixnames) <= 2)
            for i = 1:numel(ixnames)
               if ~iscell(ixnames{i})
                  ixnames{i} = {ixnames{i}};
               end
            end
         end
         % check that all index names are covered, and compute the
         % permutation
         perm = SparseTensor.get_permutation(self.components{1}.indexnames,...
                                            horzcat(ixnames{:}));
         % perm = SparseTensor.get_permutation(self.indexNames(),...
         %                                    horzcat(ixnames{:}));
         vals = self.components{1}.coefs;
         if isa(vals, 'ADI')
            % 'vals' will lose ADI status
            vals = value(vals);
         end
         
         if numel(ixnames) == 1
            ix = SparseTensor.compute_1D_index(self.components{1}.ixs(:, perm), ...
                                              shape(perm));
            M = sparse(ix, 1, vals, prod(shape), 1);
            %M = zeros(max(ix), 1);
            %M(ix) = full(self.components{1}.coefs);
         elseif numel(ixnames) == 2
            nix1 = numel(ixnames{1});
            shape1 = shape(perm(1:nix1));
            shape2 = shape(perm(nix1+1:end));
            ix1 = SparseTensor.compute_1D_index(self.components{1}.ixs(:, perm(1:nix1)), ...
                                               shape1);
            ix2 = SparseTensor.compute_1D_index(self.components{1}.ixs(:,perm(nix1+1:end)),...
                                               shape2);
            M = sparse(ix1, ix2, vals, prod(shape1), prod(shape2));
         end
         if ~force_sparse && (numel(M) < SPARSE_THRESHOLD)
            M = full(M); % convenient, for small matrices/vectors
         end
      end

      function self = expandall(self, expand_tensor)
         if nargin < 2
            expand_tensor = true; % default is true
         end
         
         if numel(self.components) == 1
            return
         end
         
         % do all contractions, start with the cheapest
         while true
            
            ixname = self.get_cheapest_pending_contraction();
            if isempty(ixname)
               % no more work to be done
               break
            end
            % do the contraction
            self = self.two_component_contraction(ixname);
            
            % if contraction led to additional contractions now being
            % associated with the same single component, carry out those too
            self = self.complete_simple_contractions();
         end
         
         % all contractions/semicontractions carried out.  Expand tensor now,
         % if requested
         if expand_tensor
            expanded_comp = self.components{1};
            for i = 2:numel(self.components)
               expanded_comp = SparseTensor.tensor_product(expanded_comp, ...
                                                          self.components{i});
               self = SparseTensor(expanded_comp);
            end
         end
      end
      
      % ---- The following methods represent implementation details, and are not ----
      % -------------------- intended to be directly run by user --------------------
      
      
      function self = complete_simple_contractions(self)
         
         for i = 1:numel(self.components)
            self.components{i} = SparseTensor.single_component_contraction(self.components{i});
         end
      end
      
      function ixname = get_cheapest_pending_contraction(self)
      
         % get all index names for which there are pending contractions
         [~, cixs] = self.indexNames();
         
         if isempty(cixs)
            ixname = [];
            return
         end
         
         % make list of all pending contractions
         cur_lowest_cost = inf;
         for cur_name = cixs
            
            cur_name = cur_name{:}; %#ok
            cost = self.contraction_cost_estimate(cur_name);
            if cost < cur_lowest_cost
               cur_lowest_cost = cost;
               ixname = cur_name;
            end
         end
      end
      

      function self = two_component_contraction(self, ixname)
         
         comp_ixs = self.component_with_ix(ixname, true);
         assert(numel(comp_ixs) == 2);
         comps = self.components(comp_ixs);
         if isa(comps{1}.coefs, 'ADI') || isa(comps{2}.coefs, 'ADI')
            self = two_component_contraction_adi(self, ixname);
         else
            self = two_component_contraction_float(self, ixname);
         end
      end
         
      
      function self = two_component_contraction_adi(self, ixname)
      
         comp_ixs = self.component_with_ix(ixname, true);
         assert(numel(comp_ixs) == 2);
         comps = self.components(comp_ixs);
         
         [c1_keep_ix, stride1, c1_contract_ix, c1_vals, c1_keep_ixnames] = ...
             SparseTensor.prepare_for_contraction(comps{1}, ixname, true);
         [c2_keep_ix, stride2, c2_contract_ix, c2_vals, c2_keep_ixnames] = ...
             SparseTensor.prepare_for_contraction(comps{2}, ixname, false);
         
         [row, col, vals] = ssparsemul([c1_keep_ix, c1_contract_ix], c1_vals,...
                                       [c2_contract_ix, c2_keep_ix], c2_vals);
         
         row = SparseTensor.compute_subs(row, stride1);
         col = SparseTensor.compute_subs(col, stride2);
          
         comp.coefs = vals;
         comp.indexnames = [c1_keep_ixnames, c2_keep_ixnames];
         comp.ixs = [row, col];
         
         keep = true(size(self.components));
         keep(comp_ixs) = false;
         self.components = [self.components(keep), comp];
      end
            
      % works well for doubles, but does not support ADI
      function self = two_component_contraction_float(self, ixname)
         
         comp_ixs = self.component_with_ix(ixname, true);
         assert(numel(comp_ixs) == 2);
         
         comps = self.components(comp_ixs);
         
         ind1 = strcmp(ixname, comps{1}.indexnames);
         rownames = comps{1}.indexnames(~ind1);
         m1 = SparseTensor(comps{1}).asMatrix({rownames,{ixname}}, true);
         m1rix = find(sum(abs(m1), 2));
         m1reduced = m1(m1rix,:);
         
         ind2 = strcmp(ixname, comps{2}.indexnames);
         colnames = comps{2}.indexnames(~ind2);
         m2 = SparseTensor(comps{2}).asMatrix({{ixname}, colnames}, true);
         m2cix = (find(sum(abs(m2), 1)))';
         m2reduced = m2(:, m2cix);
         
         %assert(issparse(m1) && issparse(m2));
         % mprod = m1 * m2;
         % size(m1reduced)
         % size(m2reduced)

         mprod = m1reduced * m2reduced;
         [i, j, v] = find(mprod);
         i = i(:); j = j(:); v = v(:);
         
         i = m1rix(i);
         j = m2cix(j);
         
         logical_size1 = max(comps{1}.ixs(:, ~ind1));
         reindex1 = cell(size(comps{1}.ixs, 2) - 1, 1);
         [reindex1{:}] = ind2sub(logical_size1, i);

         logical_size2 = max(comps{2}.ixs(:, ~ind2));
         reindex2 = cell(size(comps{2}.ixs,2)-1, 1);
         [reindex2{:}] = ind2sub(logical_size2, j);
         
         comp.coefs = v;
         comp.indexnames = [rownames, colnames];
         comp.ixs = [reindex1{:}, reindex2{:}];
         %comp.ixs = [cell2mat(reindex1'), cell2mat(reindex2')];
         
         keep = true(size(self.components));
         keep(comp_ixs) = false;
         self.components = [self.components(keep), comp];
      end
      
      
      function cost = contraction_cost_estimate(self, ixname)

         comp_ixs = self.component_with_ix(ixname, true);
         assert(numel(comp_ixs) == 2) % otherwise, contraction should already
                                      % have been carried out, since cost
                                      % estimates are only interesting for
                                      % contracting two different components
         comps = self.components;
         cost = 0;
         for i = comp_ixs
            entries = size(comps{i}.ixs, 1);
            ixind = strcmp(ixname,comps{i}.indexnames);
            numdiff = numel(unique(comps{i}.ixs(:, ixind)));
         
            cost = cost + entries / numdiff;
         end
         % comps = self.components;
         % cost = 0;
         % for i = [comp_ix_1, comp_ix_2]
         %    ixind = strcmp(ixname,comps{i}.indexnames);
         %    numdiff = numel(unique(comps{i}.ixs(:, ixind)));
            
         %    cost = cost + (size(comps{i}.ixs, 1) ^ size(comps{i}.ixs, 2) / numdiff);
         % end
      end
      
      function ixname = next_unused_contr_ix_name(self, other)
      % 'other' is an optional argument.  If provided, the produced ixname
      % should not be an existing contracting index of 'other' either.

         basename = SparseTensor.contracting_name_base();
         [~, current_contr_names] = self.indexNames();
         
         if nargin > 1
            [~, other_contr_names] = other.indexNames();
            current_contr_names = unique([current_contr_names, other_contr_names]);
         end
         
         count = 1;
         
         while true
            % search until we find an unused name
            ixname = [ basename, num2str(count)];
            if ~any(strcmp(ixname, current_contr_names))
               % we found a unique name
               return
            else
               % keep searching with a higher count value
               count = count + 1;
            end
         end
      end
      
      function comp_ix = component_with_ix(self, ixname, allow_multiple)
         if nargin < 3
            allow_multiple = false; % the usual case
         end

         res = [];
         
         for comp_ix = numel(self.components):-1:1
            if any(strcmp(ixname, self.components{comp_ix}.indexnames))
               if allow_multiple
                  res = [res, comp_ix]; %#ok
               else
                  return;
               end
            end
         end
         comp_ix = res; % multiple compoents (or zero)
      end
      
   end % end methods

   methods(Static)
      
      function [keep_ix, keep_stride, contract_ix, vals, ixnames] = ...
             prepare_for_contraction(comp, ixname, sort_contract)
         
         keep_indices = ~strcmp(ixname, comp.indexnames);
         ixnames = comp.indexnames(keep_indices);

         keep_ix = SparseTensor.compute_1D_index(comp.ixs(:, keep_indices));
         keep_stride = max(comp.ixs(:, keep_indices));
         contract_ix = comp.ixs(:,~keep_indices);

         if sort_contract
            [contract_ix, reindex] = sort(contract_ix);
            keep_ix = keep_ix(reindex);
         else
            [keep_ix, reindex] = sort(keep_ix);
            contract_ix = contract_ix(reindex);
         end
         vals = comp.coefs(reindex);
      end
      
                                                                 
      function ixname = contracting_name_base()
         ixname = 'contracting_ix__';
      end
      
      function yesno = is_contracting_ix(ixname)
         % returns true if ixname starts with the contracting name root
         cname = SparseTensor.contracting_name_base();
         yesno = (numel(cname) < numel(ixname)) && ...
                 strcmp(cname, ixname(1:numel(cname)));
      end
      
      function itensor = ind(tensor)
         itensor = tensor.toInd();
      end
         
      function comp = tensor_product(comp1, comp2)
         
         if numel(comp1.indexnames) == 0
            % comp1 is an intrinsic scalar
            comp = comp2;
            comp.coefs = comp1.coefs(1) * comp.coefs;
            return 
         elseif numel(comp2.indexnames) == 0
            % comp2 is an intrinsic scalar
            comp = comp1;
            comp.coefs = comp2.coefs(1) * comp.coefs;
            return
         else
            % both components were nontrivial tensors.  Expand tensor product
            comp.indexnames = [comp1.indexnames, comp2.indexnames];
            comp.coefs = kron(comp1.coefs, comp2.coefs);
            comp.ixs = [repelem(comp1.ixs, size(comp2.ixs, 1), 1), ...
                        repmat(comp2.ixs, size(comp1.ixs, 1), 1)];
         end
      end
      
      function comp = single_component_contraction(comp, ixname)
         if nargin < 2
            % contract all duplicated indices
            keep = cellfun(@(x) sum(strcmp(x, comp.indexnames)), ...
                           comp.indexnames) == 1;
         
            if sum(keep) == numel(comp.indexnames)
               % no index to contract
               return
            end
            comp.indexnames = comp.indexnames(keep);
            comp.ixs = comp.ixs(:, keep);
         else 
            % contract one component in one index
            local_ind = strcmp(ixname, comp.indexnames);

            % get rid of contracting index
            comp.indexnames = comp.indexnames(~local_ind);
            comp.ixs = comp.ixs(:, ~local_ind);
         end
         
         if isempty(comp.ixs)
            % result is an intrinsic scalar
            comp.coefs = sum(comp.coefs);
            comp.ixs = []; % get rid of "ghost" dimensions
            return
         end
         
         % sum up other elements at the correct place
         index1d = SparseTensor.compute_1D_index(comp.ixs);
         
         % tmp = accumarray(index1d(:), comp.coefs, [], [], [], true);
         % [uindex, ~, v] = find(tmp);
         % comp.coefs = v;
         
         % logical_size = max(comp.ixs);
         % reindex = cell(size(comp.ixs, 2), 1);
         
         % [reindex{:}] = ind2sub(logical_size, uindex);
         
         % comp.ixs = [reindex{:}];

         % keyboard
         % map = accumarray([index1d, (1:size(comp.ixs,1))'], 1, [], [], [], true);
         % keep = sum(map, 2) > 0;
         % keep_ix = find(keep);
         % elim_rows = accumarray([(1:numel(keep_ix))', keep_ix], 1, [], [], [], true);
         % map = elim_rows * map;
         % % map = map(keep,:);
         % comp.coefs = map * comp.coefs;
         % uindex = find(keep);
         
         [uindex, ~, ic] = unique(index1d);

         map = accumarray([ic, (1:size(comp.ixs, 1))'], 1, [], [], [], true);
         
         comp.coefs = map * comp.coefs;
                  
         % % extract nonzeros and recompute indices
         
         logical_size = max(comp.ixs);
         % reindex = cell(size(comp.ixs, 2), 1);
         % [reindex{:}] = ind2sub(logical_size, uindex);
         %comp.ixs = [reindex{:}];
         
         comp.ixs = SparseTensor.compute_subs(uindex, logical_size);
         % %comp.ixs = cell2mat(reindex');
      end
      
      function [t1, t2] = make_tensors_compatible(t1, t2)

         % verify that they have the same index sets
         if ~SparseTensor.is_permutation(t1.indexNames(), t2.indexNames())
            error(['Tensors cannot be made compatible as they have different ' ...
                   'index sets.']);
         end
         
         % expand both tables
         t1 = t1.expandall();
         t2 = t2.expandall();

         % make sure the order of indices is the same in both tensors
         perm = SparseTensor.get_permutation(t2.indexNames(), t1.indexNames());
         t2.components{1}.indexnames = t1.components{1}.indexnames;
         t2.components{1}.ixs = t2.components{1}.ixs(:, perm);
         
         % % fill in missing indices
         % [~, I1] = setdiff(t1.components{1}.ixs, t2.components{1}.ixs, 'rows');
         % [~, I2] = setdiff(t2.components{1}.ixs, t1.components{1}.ixs, 'rows');
         
         % t2.components{1}.ixs = [t2.components{1}.ixs; ...
         %                         t1.components{1}.ixs(I1,:)];
         % t1.components{1}.ixs = [t1.components{1}.ixs; ...
         %                         t2.components{1}.ixs(I2,:)];
         % t1.components{1}.coefs = [t1.components{1}.coefs; zeros(size(I2))];
         % t2.components{1}.coefs = [t2.components{1}.coefs; zeros(size(I1))];
         
         % % sort all indices
         % t1.components{1} = SparseTensor.sort_indices(t1.components{1});
         % t2.components{1} = SparseTensor.sort_indices(t2.components{1});
         
         % tensors should now have exactly the same indices and thus be compatible
      end
      
      function component = sort_indices(component)
         [component.ixs, I] = sortrows(component.ixs);
         component.coefs = component.coefs(I);
      end

      function tensor = apply_binary_operator(t1, t2, op)
         [t1, t2] = SparseTensor.make_tensors_compatible(t1, t2);
         assert(numel(t1.components) == 1); % should also be the case for t2 by now
         
         % determine common tensor shape
         shape = max(max(t1.components{1}.ixs), max(t2.components{1}.ixs));
         
         m1_ix1D = SparseTensor.compute_1D_index(t1.components{1}.ixs, shape);
         m2_ix1D = SparseTensor.compute_1D_index(t2.components{1}.ixs, shape);
         
         [all_ixs, ~, ic] = unique([m1_ix1D; m2_ix1D]);
         m1_ic = ic(1:numel(m1_ix1D));
         m2_ic = ic(numel(m1_ix1D)+1:end);
         n = numel(all_ixs);
         m1 = numel(m1_ix1D);
         m2 = numel(m2_ix1D);
         mat1 = accumarray([m1_ic, (1:m1)'], 1, [n, m1], [], [], true);
         mat2 = accumarray([m2_ic, (1:m2)'], 1, [n, m2], [], [], true);         
         
         v1 = mat1 * t1.components{1}.coefs;
         v2 = mat2 * t2.components{1}.coefs;
         
         res = t1.components{1};
         res.coefs = op(v1, v2);
         res.ixs = SparseTensor.compute_subs(all_ixs, shape);
         
         tensor = SparseTensor(res);
      end      
      
      % % The following function works well for doubles, but not for ADI
      % function tensor = apply_binary_operator(t1, t2, op)
      %    [t1, t2] = SparseTensor.make_tensors_compatible(t1, t2);
      %    assert(numel(t1.components) == 1); % should also be the case for t2 by now
         
      %    % determine common tensor shape
      %    shape = max(max(t1.components{1}.ixs), max(t2.components{1}.ixs));
         
      %    indices = t1.indexNames();

      %    m1 = t1.asMatrix({indices}, true, shape); % make 1D sparse vec of it
      %    m2 = t2.asMatrix({indices}, true, shape); % ditto
         
      %    if isequal(op, @rdivide)
      %       % special treatment, to avoid dividing by zeros
      %       [i, ~, v] = find(m2);
      %       m2(i) = 1./v;
      %       tmp = times(m1, m2);
      %    else
      %       tmp = op(m1, m2);
      %    end

      %    [ix, ~, val] = find(tmp);
      %    reindex = cell(numel(shape), 1);
      %    [reindex{:}] = ind2sub(shape, ix);
         
      %    res = t1.components{1};
      %    res.coefs = val(:);
      %    res.ixs = [reindex{:}];
      %    %res.ixs = cell2mat(reindex');
      %    tensor = SparseTensor(res);
      %    % res = t1.components{1};
      %    % res.coefs = op(t1.components{1}.coefs, t2.components{1}.coefs);
      %    % tensor = SparseTensor(res);
      % end
      
      function ix = compute_1D_index(multiix, shape)
         if size(multiix, 2) == 1
            ix = multiix;
            return
         end
         
         if exist('shape', 'var')
            stride = shape;
         else
            stride = max(multiix);
         end
         stride = cumprod([1, stride(1:end-1)]);
         ix = sum((multiix-1) .* stride, 2) + 1;
         
         % tmp = mat2cell(multiix, size(multiix, 1), ones(1, size(multiix, 2)));
         % ix = sub2ind(max(multiix), tmp{:});
      end
      
      function ixs = compute_subs(ix1d, shape)
         
         stride = cumprod([1, shape(1:end-1)]);
                  
         tmp = bsxfun(@(x,y) ceil(x/y), ix1d, stride);
         ixs = bsxfun(@mod, tmp-1, shape)+1;
      end
      
      
      function isperm = is_permutation(cells1, cells2)
         isperm = ...
             (numel(cells1) == numel(cells2)) && ...
             isempty(setdiff(cells1, cells2)) && ...
             isempty(setdiff(cells2, cells1));
      end
      
      function perm = get_permutation(cells1, cells2)
         % check that this actually is a permutation
         if ~SparseTensor.is_permutation(cells1, cells2)
            error('provided indices do not constitute a permutation of tensor indices.');
         end
         perm = zeros(1, numel(cells1));
         for ix = 1:numel(cells1)
            perm(strcmp(cells1{ix}, cells2)) = ix;
         end
      end

      function components = make_matrix_tensor(coefs, indexnames)

         % if user did not wrap the index name in a cell, do it here
         if ~iscell(indexnames)
            indexnames = {indexnames};
         end
         
         if isscalar(coefs) 
            component.coefs = coefs;
            component.ixs = [];
            assert(numel(indexnames) == 0);
         elseif isvector(coefs)
            nz = find(coefs);
            component.coefs = coefs(nz);
            component.coefs = component.coefs(:); % ensure column vec
            component.ixs = nz(:);
            assert(numel(indexnames) == 1);
         else
            assert(ismatrix(coefs));
            if issparse(coefs)
               nz = find(coefs); % only keep nonzeros
            else
               nz = (1:numel(coefs))';
            end
            [i, j] = ind2sub(size(coefs), nz);
            component.coefs = coefs(nz);
            component.coefs = component.coefs(:); % ensure column vec
            component.ixs = [i, j];
         end
         component.indexnames = indexnames(:)';         
         components = {component};
      end
   end
end % end classdef
