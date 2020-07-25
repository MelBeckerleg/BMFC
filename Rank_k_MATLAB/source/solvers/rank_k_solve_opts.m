classdef rank_k_solve_opts

properties
    it %use AM to update rank1?
    rank_1_method %for a recovery method
    recovery_solver %for an impute and cluster method
    k %rank to solve for
    pf %in case of post-processing (columnwise partitioning) step 
        %note that this does not change the observed behaviour
    %for Neighbour/Nuclear
    binarised
    recovered
    epsilon %for 1bitMC
    thresh %for 1bitMC
end
end
