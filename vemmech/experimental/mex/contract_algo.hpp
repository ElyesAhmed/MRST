#ifndef __CONTRACT_ALGO_HPP
#define __CONTRACT_ALGO_HPP

#include <set>
#include <numeric>
#include <algorithm>
#include <array>
#include <assert.h>

#include "TensorComp.hpp"

// ----------------------------------------------------------------------------
template<typename T> inline std::vector<std::string>
identify_contracting_multiindex(const std::vector<TensorComp<T>>& comps)
// ----------------------------------------------------------------------------
{
  // extract all indices
  std::vector<std::string> all_ixs;
  for (const auto& c : comps)
    std::copy(c.indexNames().begin(),
              c.indexNames().end(),
              std::back_inserter(all_ixs));

  // identify repeated indices
  std::set<std::string> contr_ixs;
  std::sort(all_ixs.begin(), all_ixs.end());
  for (int i = 0; i != all_ixs.size()-1; ++i)
    if (all_ixs[i] == all_ixs[i+1])
      contr_ixs.insert(all_ixs[i]);

  // converting into a vector
  const std::vector<std::string> cixs(contr_ixs.begin(), contr_ixs.end());

  
  // determine order of indices (start index should have highest number of
  // unique values)
  std::vector<size_t> ix_numvals(cixs.size(), 0);
  for (int i = 0; i != cixs.size(); ++i)
    for (const auto& c : comps)
      ix_numvals[i] = std::max(ix_numvals[i], c.numUniqueValuesFor(cixs[i]));

  std::vector<std::array<int, 2>> sort_ix;
  for (int i = 0; i != ix_numvals.size(); ++i)
    sort_ix.emplace_back( std::array<int, 2>{(int)ix_numvals[i], i} );

  std::sort(sort_ix.begin(), sort_ix.end());
  std::reverse(sort_ix.begin(), sort_ix.end());

  std::vector<int> perm(sort_ix.size());
  std::transform(sort_ix.begin(), sort_ix.end(), perm.begin(),
                 [](const std::array<int, 2>& el) {return el[1];});
    
  std::vector<std::string> sorted_ixs(cixs.size());
  for (int i = 0; i != cixs.size(); ++i)
    sorted_ixs[i] = cixs[perm[i]];
  //return cixs;
  return sorted_ixs;
  
}

// ----------------------------------------------------------------------------
template<typename T> std::vector<int>
identify_comps_for_multiix_iteration(const std::vector<TensorComp<T>>& comps,
                                     const std::vector<std::string>& cixs)
// ----------------------------------------------------------------------------
{
  std::vector<int> result(cixs.size(), -1); // -1 flags that index not yet set
  std::vector<int> count(cixs.size(), -1);

  for (int c_ix = 0; c_ix != comps.size(); ++c_ix) {
    for (int i = 0; i != cixs.size(); ++i) {
      const size_t numvals = comps[c_ix].numUniqueValuesFor(cixs[i]);
      if (numvals > 0 && (count[i] < 0 || count[i] > numvals)) {
        result[i] = c_ix;
        count[i] = numvals;
      }
    }
  }
  return result;
}

// ----------------------------------------------------------------------------
inline void advance_multiindex(std::vector<size_t>& mix, const std::vector<size_t>& bounds)
// ----------------------------------------------------------------------------
{
  mix.back() = mix.back() + 1;
  for (int i = mix.size() - 1; i > 0; --i) {
    if (mix[i] == bounds[i]) {
      mix[i] = 0;
      mix[i-1] = mix[i-1] + 1;
    }
  }
  // std::cout << "new mix: ";
  // for (int i = 0; i != mix.size(); ++i)
  //   std::cout << mix[i] << "  (" << bounds[i] << ") ";
  // std::cout <<  std::endl;
}

// ----------------------------------------------------------------------------
template<typename Index>
std::vector<size_t> sort_multiindex(std::vector<Index>& mix, int num_elem)
// ----------------------------------------------------------------------------
{
  struct SortElem {
    const Index* ptr;
    size_t init_pos;
  };
  size_t N = mix.size() / num_elem;
  std::vector<SortElem> mixptrs(N);
  for (size_t i = 0; i != N; ++i)
    mixptrs[i] = { &mix[i * num_elem], i};

  auto comp = [num_elem](const SortElem& i1, const SortElem& i2) {
                for (int i = 0; i != num_elem; ++i)
                  if (i1.ptr[i] != i2.ptr[i])
                    return i1.ptr[i] < i2.ptr[i];
                return false;
              };
  std::sort(mixptrs.begin(), mixptrs.end(), comp);

  // copying sorted elements back
  std::vector<size_t> reindexes(N);
  std::vector<Index> result(mix.size());
  size_t count = 0;
  size_t reindex_count = 0;
  for (auto it = mixptrs.begin(); it != mixptrs.end(); ++it) {
    for (int i = 0; i != num_elem; ++i) 
      result[count++] = it->ptr[i];
    
    reindexes[reindex_count++] = it->init_pos;
  }
  std::swap(mix, result);
  return reindexes;
}

// ----------------------------------------------------------------------------
template<typename Index>
void transpose(std::vector<Index>& vals, size_t num_cols)
// ----------------------------------------------------------------------------
{
  const size_t num_rows = vals.size() / num_cols;
  std::vector<Index> result(vals.size());
  size_t ix = 0;
  for(size_t i = 0; i != num_rows; ++i)
    for(size_t j = 0; j != num_cols; ++j)
      result[ix++] = vals[j * num_rows + i];
  std::swap(vals, result);
}

// ----------------------------------------------------------------------------
template<typename T>
std::vector<typename TensorComp<T>::Index>
compute_multiindices(
     const std::vector<std::string>& cixs,
     const std::vector<int>& mix_comp,
     const std::vector<TensorComp<T>>& comps)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  typedef std::array<Index, 2> RangeEntry;

  struct IxControl {
    int comp_ix;
    std::vector<Index> ix_entries;
    std::vector<std::vector<Index>> ix_values;
  };
  std::vector<IxControl> controls;
  for (int c_ix = 0; c_ix != comps.size(); ++c_ix) {

    // determine indices for which this component will control the multiindex
    std::vector<std::string> ixnames;
    std::vector<Index> ix_pos;
    for (int i = 0; i != mix_comp.size(); ++i) {
      if (mix_comp[i] == c_ix) {
        ixnames.push_back(cixs[i]);
        ix_pos.push_back(i);
      }
    }
    if (ixnames.empty())
      continue;
    
    std::vector<Index> ixvals;
    for (int i = 0; i != ixnames.size(); ++i) {
      const auto col = comps[c_ix].indexValuesFor(ixnames[i]);
      ixvals.insert(ixvals.end(), col.begin(), col.end());
    }
    transpose(ixvals, ixnames.size());
    sort_multiindex(ixvals, ixnames.size());
    
    // identify where index changes occur
    const size_t elnum = ixnames.size();
    const size_t N = ixvals.size() / elnum;
    
    std::vector<Index> changes(1, 0);
    for (int i = 1; i != N; ++i) {
      for (int j = 0; j != elnum && changes.back() != i; ++j)
        if (ixvals[elnum * i + j] != ixvals[elnum * (i-1) + j])
          changes.push_back(i);
    }

    std::vector<std::vector<Index>> unique_ixvals;
    for (int i = 0; i != elnum; ++i) {
      std::vector<Index> unique_ival(changes.size());

      for (int j = 0; j != changes.size(); ++j)
        unique_ival[j] = ixvals[changes[j] * elnum + i];

      unique_ixvals.push_back(unique_ival);
    }
    controls.emplace_back(IxControl {c_ix, ix_pos, unique_ixvals} );
  }

  // sort controls so that highest stride lies in first index
  std::sort(controls.begin(), controls.end(),
            [](const IxControl& c1, const IxControl& c2) {
              return c1.ix_entries[0] < c2.ix_entries[0];});
  
  std::vector<size_t> indices(controls.size(), 0); 
  std::vector<size_t> indices_max;
  for (const auto& c : controls)
    indices_max.push_back(c.ix_values[0].size());
  assert(indices.size() > 0);
  
  std::vector<Index> result;
  std::vector<Index> new_multiindex(mix_comp.size());
  
  while (indices[0] != indices_max[0]) {
    for (int i = 0; i != controls.size(); ++i) {
      const auto& ctrl = controls[i];
      for (int j = 0; j != ctrl.ix_values.size(); ++j)
        new_multiindex[ctrl.ix_entries[j]] = ctrl.ix_values[j][indices[i]];
    }
    result.insert(result.end(), new_multiindex.begin(), new_multiindex.end());
    
    // advance loop multiindex
    advance_multiindex(indices, indices_max);
  }
  sort_multiindex(result, mix_comp.size());
  
  return result;
}

// ----------------------------------------------------------------------------
template<typename T> std::vector<std::array<typename TensorComp<T>::Index, 2>>
compute_index_ranges(const TensorComp<T>& comp,
                     const std::vector<std::string>& cixnames,
                     const std::vector<typename TensorComp<T>::Index>& mixs)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  const int elnum = (int)cixnames.size(); // number of multiindex components
  const size_t mixnum = mixs.size() / elnum; // number of multiindices
  
  // determine contracting indices present in this component
  std::vector<int> active;
  for (int i = 0; i != elnum; ++i)
    if (comp.indexPos(cixnames[i]) >= 0)
      active.push_back(i);

  const int actnum = active.size(); // number of active contracting indices in comp

  // preparing index vectors for comparison
  std::vector<std::vector<Index>> tmp;
  for (size_t i = 0; i != actnum; ++i)
    tmp.push_back(comp.indexValuesFor(cixnames[active[i]]));

  std::vector<Index> comp_ixs;
  for (size_t i = 0; i != tmp[0].size(); ++i) 
    for (int j = 0; j != tmp.size(); ++j)
      comp_ixs.push_back(tmp[j][i]);

  std::vector<Index> mix_ixs;
  for (size_t i = 0; i != mixnum; ++i)
    for (int j = 0; j != actnum; ++j)
      mix_ixs.push_back(mixs[i*elnum + active[j]]);

  const auto reindex = sort_multiindex(mix_ixs, actnum);
                                       

  // loop through all multiindices and define ranges
  std::vector<std::array<Index, 2>> result(mixnum);

  const auto ixless = [actnum](typename std::vector<Index>::iterator cix,
                               typename std::vector<Index>::iterator mix) {
    for (int i = 0; i != actnum; ++i)
      if (cix[i] != mix[i])
        return (cix[i] < mix[i]);
    return false;
  };

  auto cit = comp_ixs.begin();

  for (auto mit = mix_ixs.begin(); mit != mix_ixs.end(); mit += actnum) {
    
    const int ix = (mit - mix_ixs.begin()) / actnum;

    // check if this is a repeat of last iterations values for mit
    if (ix > 0 && !( ixless(mit, mit-actnum) || ixless(mit-actnum, mit))) {
      result[ix] = result[ix-1];
      continue;
    }
          
    // advance counter until multiindex is found
    while (ixless(cit, mit) && cit < comp_ixs.end())
      cit += actnum;

    if (ixless(mit, cit)) // indices are not equal
      continue;

    // if we got here, the two indices are equal
    result[ix][0] = (cit - comp_ixs.begin()) / actnum;

    while (!ixless(mit, cit) && cit < comp_ixs.end())
      cit += actnum;

    result[ix][1] = (cit - comp_ixs.begin()) / actnum;

  }
  // reindex results according to original order of elements in mixs
  std::vector<std::array<Index, 2>> result_reorder(result.size());
  for (size_t i = 0; i != result.size(); ++i)
    result_reorder[reindex[i]] = result[i];

  std::swap(result, result_reorder);
  
  return result;
}

// ----------------------------------------------------------------------------
template<typename T>
std::pair<std::vector<typename TensorComp<T>::Index>, std::vector<std::string>>
free_indices_for(const TensorComp<T>& comp,
                 const std::vector<std::string>& cixnames)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;

  const auto ixnames = comp.indexNames();
  std::vector<std::pair<std::vector<Index>, std::string>> free_ixs;

  for (const auto& iname : ixnames)
    if (std::find(cixnames.begin(), cixnames.end(), iname) == cixnames.end())
      // this is a free index, since it is not among the contracting indices
      free_ixs.emplace_back( std::pair<std::vector<Index>, std::string> {
          comp.indexValuesFor(iname), iname
            });;

  if (free_ixs.empty())
    return std::pair<std::vector<Index>, std::vector<std::string>>();
  
  // reshuffle indices
  std::vector<Index> res_ixs;
  const size_t numvals = free_ixs[0].first.size();
  for (Index i = 0; i != numvals; ++i)
    for (int j = 0; j != free_ixs.size(); ++j)
      res_ixs.push_back(free_ixs[j].first[i]);

  // collecting free index names
  std::vector<std::string> free_ix_names;
  for (const auto& fix : free_ixs)
    free_ix_names.push_back(fix.second);

  return std::pair<std::vector<Index>, std::vector<std::string>> {
    res_ixs, free_ix_names
  };
}

// ----------------------------------------------------------------------------
template<typename T> TensorComp<T>
compute_sums(const std::vector<TensorComp<T>>& comps,
             const std::vector<std::vector<std::array<typename TensorComp<T>::Index, 2>>>& ranges,
             const std::vector<std::string>& cixnames)
// ----------------------------------------------------------------------------
{
  typedef typename TensorComp<T>::Index Index;
  // setup free indices associated with each component
  std::vector<std::pair<std::vector<Index>, std::vector<std::string>>> free_indices;
  std::vector<int> free_indices_num;
  for (const auto& c : comps) {
    free_indices.emplace_back(free_indices_for(c, cixnames));
    free_indices_num.push_back(free_indices.back().second.size());
  }
  const size_t N = ranges[1].size();
  std::vector<T> coefs;
  std::vector<Index> indices;
  std::vector<size_t> rstart(ranges.size()), rlen(ranges.size()), running(ranges.size(), 0);
  
  for (size_t i = 0; i != N; ++i) {

    bool skip = false;
    for (int j = 0; j < ranges.size() && !skip; ++j) {
      running[j] = 0;
      rstart[j] = ranges[j][i][0];
      rlen[j] = ranges[j][i][1] - rstart[j];
      skip = (rlen[j] == 0);
    }
    if (skip)
      continue; // no (value) for this element
    
    while (running[0] != rlen[0]) {

      T new_coef = 1;

      for (int j = 0; j != ranges.size(); ++j) {
        const Index cur_ix = rstart[j] + running[j];
        const Index cur_ix2 = cur_ix * free_indices_num[j];
        
        new_coef *= comps[j].coefs()[cur_ix];
        
        indices.insert(indices.end(),
                       &(free_indices[j].first[cur_ix2]),
                       &(free_indices[j].first[cur_ix2 + free_indices_num[j]]));
      }
      
      coefs.push_back(new_coef);
      
      advance_multiindex(running, rlen);
    }

  }

  // collect the index names of the to-be-created tensor
  std::vector<std::string> indexnames;
  for (const auto& fi : free_indices)
    indexnames.insert(indexnames.end(), fi.second.begin(), fi.second.end());

  // reformat indices (unless tensor is an intrinsic scalar, in which case we return it)
  const size_t num_free_indices = indexnames.size();
  if (num_free_indices == 0)
    // tensor is an intrinsic scalar
    return TensorComp<T>(std::accumulate(coefs.begin(), coefs.end(), (T)0));

  const size_t num_free_ix_values = indices.size() / num_free_indices;

  
  std::vector<Index> indices_reordered(indices.size());
  size_t pos = 0;
  for (size_t j = 0; j != num_free_indices; ++j)
    for (size_t i = 0; i != num_free_ix_values; ++i)
      indices_reordered[pos++] = indices[i * num_free_indices + j];
  
  // generate result tensor and add up element with similar indices
  TensorComp<T> result(indexnames, coefs, indices_reordered);
  result.sortIndicesByNumber(true).sortElementsByIndex().sumEqualIndices();

  return result;

}

// ============================================================================
template<typename T>
std::vector<TensorComp<T>> contract_components(std::vector<TensorComp<T>> comps)
// ============================================================================
{
  using std::cout;
  using std::endl;
  typedef typename TensorComp<T>::Index Index;
  // identify contracting indices; sort them according to number of unique values
    const std::vector<std::string> cixs(identify_contracting_multiindex(comps));

  // Determine which components will control the iteration of each index in the
  // multiindex.  
  const auto mix_comp = identify_comps_for_multiix_iteration(comps, cixs);

  // prepare components for iteration  by sorting index order and coefficient properly
  for (int i = 0; i != comps.size(); ++i)  
    comps[i].sortIndicesByNumber(true).moveIndicesFirst(cixs).sortElementsByIndex();
    
  // identify all possible multiindex values and corresponding entry ranges in
  // components
  const std::vector<Index> mix = compute_multiindices(cixs, mix_comp, comps);

  std::vector<std::vector<std::array<Index, 2>>> comp_ranges;

  for (const auto& c : comps)
    comp_ranges.emplace_back(compute_index_ranges(c, cixs, mix));

  // multiply together relevant indices and construct new component
  TensorComp<T> summed_tensor = compute_sums(comps, comp_ranges, cixs);

  std::vector<TensorComp<T>> result(1, summed_tensor);
  return result;
  //return comps; // @@ implement properly
}
  
#endif
