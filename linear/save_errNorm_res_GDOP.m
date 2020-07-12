function [ned_err_solver,hor_err_solver,err_solver,GDOP_solver,nsv,dnsv,nprior,dnprior,res_solver] =...
    save_errNorm_res_GDOP...
    (solver_n,Grdpose,re_pos_solver,GDOP,nsv,dnsv,nprior,dnprior,max_num_sv,ind_mark,solverRes_cells)

    err_solver     = NaN(solver_n,1); 
    hor_err_solver = NaN(solver_n,1);    
    ned_err_solver = NaN(solver_n,1);
    GDOP_solver    = NaN(solver_n,1);
    res_solver     = [];
    
    for idx=1:solver_n
        [pos_llh,~,~]=ecef2llh_iter(re_pos_solver{idx});
        R_e2g =ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
        err_pos = Grdpose - re_pos_solver{idx}; % position err = true - estimated position
        ned_err=R_e2g*err_pos;    
        ned_err_solver(idx) = norm(ned_err);
        hor_err_solver(idx) = norm(ned_err(1:2));
        err_solver(idx)     = norm(err_pos);
        GDOP_solver(idx)    = GDOP(idx,1);
        %------------------------------%
        tmp_res = NaN(max_num_sv,1);
        tmp_res(ind_mark) = solverRes_cells{idx};
        res_solver = [res_solver; tmp_res];
    end
end

