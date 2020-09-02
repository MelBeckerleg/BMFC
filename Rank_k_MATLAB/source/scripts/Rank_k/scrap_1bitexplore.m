addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/error'))
addpath(genpath('/home/user/mosek/'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/solvers/MatrixCompletion'))
clear all
tag='aug_1bitmed_run6';
tag_test_params = 'num_methods=3;m_vals=[500];k_vals=[0.1:0.05:0.3];rho_vals=[0.1:0.1:0.7];phase_type=1;num_trials=20;tau=0.3;epsilon=0.05;';%this is for my memory, nothing is set from this!
eval(tag_test_params);m=m_vals(1);n=m;
synthetic_set_up='Brickwork';
%rho_vals = [0.1 0.2 0.3]; 
%k_vals = [0.1 0.2 0.3];
names={{'1bit'},{'1bitR'},{'1bitRT'}}
RESULTS.err_hamming = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.train_err = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.test_err = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));

RESULTS.row_cluster_strength = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.col_cluster_strength = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.recovered_rank = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.num_rows_recovered = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.num_rows_recovered_95 = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.num_rows_recovered_50 = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.num_rows_recovered_80 = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));
RESULTS.num_rows_recovered_75 = zeros(num_trials,num_methods,length(k_vals),length(rho_vals));



rho_idx = 0; k_idx =0 ; 
for rho = rho_vals
    rho_idx=rho_idx+1
    exact_rho = rho;
    for k_val = k_vals 
        k_idx = k_idx+1
        k = ceil(m*k_val);

for trial_idx=[1:num_trials]

rng(trial_idx+50);
epsilon=0.05;
prob=rank_k_problem();
prob.m=m;prob.n=n;prob.epsilon=epsilon;prob.rho=rho;prob.tau=tau;prob.k=k;prob.k_solve=k;
prob.set_up=synthetic_set_up; % generate_row_clusters
[prob.X_true,prob.Y_true,prob.A_true,prob.W_omega]= generate_row_clusters(prob);
prob.mask=1.*(abs(prob.W_omega)>0);

opts.epsilon= epsilon;

        if 0

            opts.k=k;opts.recovered=0;opts.binarised=1;
            opts.recovery_solver='1bitMC';
            [Xr,Yr] = Recover(prob.W_omega,opts,prob.A_true);
        end

    if 1

        CompletedMat = OneBitMC(m,n,k,[],prob.W_omega,opts);
        %set the threshold according to the definition in chapter 5,
        c_rho = n^(1/24)*exact_rho;
        c_k = n^(1/12)*k;
        %take sqrt(T_epsilon) since pdist calculates Euclidean
        %distance:
        opts.thresh = sqrt(32*sqrt(2)*sqrt(c_k/c_rho)*(1+1/(4*(1/2-opts.epsilon)^2)) ...
                                   *(1/2-opts.epsilon)*(1+(1/2-opts.epsilon)^2/(1/4-(1/2-opts.epsilon)^2)) ...
                                   *8*sqrt(2)*(1+sqrt(6))/exp(1) ...
                                   *sqrt(1+(m+n)*log(m)/(exact_rho*m*n)) ...
                                   *n^(3/4-1/6+1/24)); %this is far too large!
        opts.thresh=sqrt(2*sqrt(opts.thresh)); %for smaller m, use this as a proxy
        %opts.thresh=sqrt(n^(3/4-1/6+1/24));

        sprintf('The error is \n')
        val = nnz(CompletedMat - prob.A_true)/(n*m)


        sprintf('The binarised error is \n')
        val = nnz((CompletedMat>0.0)*1. - prob.A_true)/(n*m)

        recover_rows_function = @threshold_group;%don't bother with nearest, 
                                                 %as the recovery region 
                                                 %matches so well! recover_rows_function = @nearest_group;
        X = recover_rows_function((CompletedMat>0)*1.,k,opts);
        %before you next run: set up to do thi
        Yr1 = majority_vote(CompletedMat,X);
        Yr2 = majority_vote((CompletedMat>0)*2.-1,X);
        Yobs = majority_vote(prob.W_omega,X);

    end

trues = prob.A_true; 


method_idx = 0;
for Y_tag = [{'Yobs'},{'Yr1'}, {'Yr2'}]
    method_idx= method_idx+1;
    Y = eval(Y_tag{1});    
    preds = (X*Y'>0)*1.;
    sprintf('The final error is \n')
    RESULTS.err_hamming(trial_idx,method_idx,k_idx,rho_idx) = nnz(preds - prob.A_true)/500^2
    RESULTS.err_train(trial_idx,method_idx,k_idx,rho_idx) = nnz(preds(find(prob.W_omega)) - trues(find(prob.W_omega)))/nnz(prob.W_omega);
    RESULTS.err_test(trial_idx,method_idx,k_idx,rho_idx) = nnz(preds(find(1 - abs(prob.W_omega))) - trues(find(1- abs(prob.W_omega))))/(500^2-nnz(prob.W_omega));
    [pcerr,pcerr_col] = per_cluster_error(X,Y,prob.A_true) ; 
    RESULTS.row_cluster_strength(trial_idx,method_idx,k_idx,rho_idx) = pcerr;
    RESULTS.col_cluster_strength(trial_idx,method_idx,k_idx,rho_idx) = pcerr_col;
    RESULTS.recovered_rank(trial_idx,method_idx,k_idx,rho_idx) = min(length(find(sum(Y))),length(find(sum(X))));
    RESULTS.num_rows_recovered(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0);
    RESULTS.num_rows_recovered_95(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.05);
    RESULTS.num_rows_recovered_50(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.5);
    RESULTS.num_rows_recovered_80(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.2);
    RESULTS.num_rows_recovered_75(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.25);
end


fn=sprintf('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_%s_tau_%s_eps_%s_%s.mat',string(m),string(100*tau),string(100*epsilon),tag);
save(fn,'names','tag','RESULTS','tag_test_params');





end
    end
    k_idx = 0;
end


function [proportion] = row_recovered(A,B,thresh)
    [m,n] = size(A);
    prop_different_by_row = 1/n*sum(abs(A-B),2);
    proportion = 1/m*sum(prop_different_by_row <= thresh);
end