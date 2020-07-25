classdef rank_k_problem
    properties
        set_up='generated';
        m    %number of rows
        n    %number of columns
        k    %rank
        tau  %proportion of positives
        tau1 %proportion of positives by row (used for rank 1 planted tile model)
        tau2 %proportion of positives by column (used for rank 1 planted tile model)
        rho     %subampling factor
        epsilon %noise
        A_true  %Ground truth, binary matrix
        mask    %mask of observed entries 
        W_omega  %Noisy observed matrix, +1=postive, -1:= negative and 0:= unobserved
        X_true    %ground row cluster
        Y_true    %ground truth column cluster
        X    %approx row cluster
        Y    %approx column cluster
    end
    
    methods
        generate_function %generative model
        load_function    %load in existing dataset   
        %nosie/subsampling is contained in the above two - for now... 
    end
    
    properties
        k_solve %rank of the approximation
    end
end

