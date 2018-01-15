//
// include necessary system headers
//
#include <mex.h>
#include "matrix.h"

#include <iostream>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/runtime.hpp>


#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/zero_copy.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/runtime.hpp>

#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/preconditioner/runtime.hpp>




#include <string>


#include <amgcl/preconditioner/cpr.hpp>
#include <amgcl/preconditioner/cpr_drs.hpp>

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *result; 
    double *rhs;
    double *err;
    mwSize m,n,nnz;
    mwIndex * cols;
    mwIndex * rows;
    const mxArray * pa;

    double * entries;
    std::string relaxParam;
    std::string coarsenParam;
    
    if (nrhs != 3) { 
	    mexErrMsgTxt("3 input arguments required."); 
    } else if (nlhs > 2) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 

    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||  !mxIsSparse(prhs[0]) ) { 
	    mexErrMsgTxt("Matrix should be a real sparse matrix."); 
        return;
    } 
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) { 
	    mexErrMsgTxt("Right hand side must be real double column vector.");
        return; 
    } 
    // main();
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    result = mxGetPr(plhs[0]);
    err = mxGetPr(plhs[1]);
    
    cols    = mxGetJc(prhs[0]);
    rows    = mxGetIr(prhs[0]);
    entries = mxGetPr(prhs[0]);
    nnz  = mxGetNzmax(prhs[0]);
    rhs     = mxGetPr(prhs[1]);
    pa = prhs[2];
    double tolerance = mxGetScalar(mxGetField(pa, 0, "tolerance"));
    int maxiter = mxGetScalar(mxGetField(pa, 0, "maxIterations"));
    int block_size = mxGetScalar(mxGetField(pa, 0, "block_size"));
    bool use_drs = mxGetScalar(mxGetField(pa, 0, "use_drs"));
    
    
    int coarsen_id = mxGetScalar(mxGetField(pa, 0, "coarsening"));
    int relax_p_id = mxGetScalar(mxGetField(pa, 0, "relaxation"));
    int relax_s_id = mxGetScalar(mxGetField(pa, 0, "s_relaxation"));
    int solver_id = mxGetScalar(mxGetField(pa, 0, "solver"));
    
    // int precond_id = mxGetScalar(mxGetField(pa, 0, "preconditioner"));
    
    std::vector<double> b(n);
    for(int ix = 0; ix < n; ix++){
        b[ix] = rhs[ix];
    }
    int M = (int)m;
    /***************************************
     * Start AMGCL-link and select options *
     ***************************************/
    typedef amgcl::backend::builtin<double> Backend;
    
    typedef
        amgcl::runtime::amg<Backend>
        PPrecond;

    typedef
        amgcl::runtime::relaxation::as_preconditioner<Backend>
        SPrecond;


    
    boost::property_tree::ptree prm;
    /* Set tolerance */
    prm.put("solver.tol", tolerance);
    if(maxiter > 0){
        prm.put("solver.maxiter", maxiter);
    }
    prm.put("precond.block_size", block_size);

    /* Select coarsening strategy */
    coarsenParam = "precond.pprecond.coarsening.type";
    switch(coarsen_id) {
        case 1: 
            prm.put(coarsenParam,  amgcl::runtime::coarsening::smoothed_aggregation);
            break;
        case 2: 
            prm.put(coarsenParam,  amgcl::runtime::coarsening::ruge_stuben);
            break;
        case 3: 
            prm.put(coarsenParam,  amgcl::runtime::coarsening::aggregation);
            break;
        case 4: 
            prm.put(coarsenParam,  amgcl::runtime::coarsening::smoothed_aggr_emin);
            break;
        default : mexErrMsgTxt("Unknown coarsen_id."); 
    }
    /* When is a level coarse enough */
    int coarse_enough = mxGetScalar(mxGetField(pa, 0, "coarse_enough"));
    if (coarse_enough >= 0){
        prm.put("precond.pprecond.coarse_enough", coarse_enough);
    }
    /* Use direct solver for coarse sys */
    bool direct_coarse = mxGetScalar(mxGetField(pa, 0, "direct_coarse"));
    prm.put("precond.pprecond.direct_coarse", direct_coarse);
    /* Max levels */
    int max_levels = mxGetScalar(mxGetField(pa, 0, "max_levels"));
    if (max_levels >= 0){
        prm.put("precond.pprecond.max_levels", max_levels);
    }
    /* Number of cycles */
    int ncycle = mxGetScalar(mxGetField(pa, 0, "ncycle"));
    if (ncycle >= 0){
        prm.put("precond.pprecond.ncycle", ncycle);
    }
    /* Pre cycles */
    int npre = mxGetScalar(mxGetField(pa, 0, "npre"));
    if (npre >= 0){
        prm.put("precond.pprecond.npre", npre);
    }
    /* Post cycles */
    int npost = mxGetScalar(mxGetField(pa, 0, "npost"));
    if (npost >= 0){
        prm.put("precond.pprecond.npost", npost);
    }
    /* Pre cycles (precond) */
    int pre_cycles = mxGetScalar(mxGetField(pa, 0, "pre_cycles"));
    if (pre_cycles >= 0){
        prm.put("precond.pprecond.pre_cycles", pre_cycles);
    }

    /* Select relaxation strategy */
    relaxParam = "precond.pprecond.relax.type";
    switch(relax_p_id) {
        case 1: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::spai0);
            break;
        case 2: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::gauss_seidel);
            break;
        case 3: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::ilu0);
            break;
        case 4: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::iluk);
            break;
        case 5: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::ilut);
            break;
        case 6: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::damped_jacobi);
            break;
        case 7: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::spai1);
            break;
        case 8: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::chebyshev);
            break;
        default : mexErrMsgTxt("Unknown relax_id."); 
    }
    
    relaxParam = "precond.sprecond.type";
    switch(relax_s_id) {
        case 1: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::spai0);
            break;
        case 2: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::gauss_seidel);
            break;
        case 3: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::ilu0);
            break;
        case 4: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::iluk);
            break;
        case 5: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::ilut);
            break;
        case 6: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::damped_jacobi);
            break;
        case 7: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::spai1);
            break;
        case 8: 
            prm.put(relaxParam,  amgcl::runtime::relaxation::chebyshev);
            break;
        default : mexErrMsgTxt("Unknown relax_id."); 
    }
    /* Select solver */
    switch(solver_id) {
        case 1: 
            prm.put("solver.type",  amgcl::runtime::solver::bicgstab);
            break;
        case 2: 
            prm.put("solver.type",  amgcl::runtime::solver::cg);
            break;
        case 3: 
            prm.put("solver.type",  amgcl::runtime::solver::bicgstabl);
            {
                int L = mxGetScalar(mxGetField(pa, 0, "bicgstabl_l"));
                prm.put("solver.L", L);
                double delta = mxGetScalar(mxGetField(pa, 0, "bicgstabl_delta"));
                prm.put("solver.delta", delta);
                bool convex = mxGetScalar(mxGetField(pa, 0, "bicgstabl_convex"));
                prm.put("solver.convex", convex);
            }
            break;
        case 4: 
            prm.put("solver.type",  amgcl::runtime::solver::gmres);
            {
                int M = mxGetScalar(mxGetField(pa, 0, "gmres_m"));
                prm.put("solver.M", M);
            }
            break;
        case 5: 
            prm.put("solver.type",  amgcl::runtime::solver::lgmres);
            {
                int M = mxGetScalar(mxGetField(pa, 0, "gmres_m"));
                prm.put("solver.M", M);
                int K = mxGetScalar(mxGetField(pa, 0, "lgmres_k"));
                prm.put("solver.K", K);
                bool always_reset = mxGetScalar(mxGetField(pa, 0, "lgmres_always_reset"));
                prm.put("solver.always_reset", always_reset);
                bool store_Av = mxGetScalar(mxGetField(pa, 0, "lgmres_store_av"));
                prm.put("solver.store_Av", store_Av);            
            }
            break;
        case 6: 
            prm.put("solver.type",  amgcl::runtime::solver::fgmres);
            {
                int s = mxGetScalar(mxGetField(pa, 0, "idrs_s"));
                prm.put("solver.s", s);
                double omega = mxGetScalar(mxGetField(pa, 0, "idrs_omega"));
                prm.put("solver.omega", omega);
                bool replace = mxGetScalar(mxGetField(pa, 0, "idrs_replacement"));
                prm.put("solver.replacement", replace);
            }
            break;
        case 7: 
            prm.put("solver.type",  amgcl::runtime::solver::idrs);
            break;
        default : mexErrMsgTxt("Unknown solver_id."); 
    }


    /***************************************
     *        Solve problem                *
     ***************************************/
    std::vector<double> x(M, 0.0);
    int    iters;
    double error;
    
    if(use_drs){
        double dd = mxGetScalar(mxGetField(pa, 0, "drs_eps_dd"));
        double ps = mxGetScalar(mxGetField(pa, 0, "drs_eps_ps"));
        prm.put("precond.eps_dd", dd);
        prm.put("precond.eps_ps", dd);
        amgcl::make_solver<
            amgcl::preconditioner::cpr_drs<PPrecond, SPrecond>,
            amgcl::runtime::iterative_solver<Backend>
            > solve(amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]), prm);
            boost::tie(iters, error) = solve(b, x);
    }else{
        amgcl::make_solver<
            amgcl::preconditioner::cpr<PPrecond, SPrecond>,
            amgcl::runtime::iterative_solver<Backend>
            > solve(amgcl::adapter::zero_copy(n, &cols[0], &rows[0], &entries[0]), prm);
            boost::tie(iters, error) = solve(b, x);
    }
    for(int ix=0; ix < M; ix++){
        result[ix] = x[ix];
    }
    err[0] = error;
    return;
}

