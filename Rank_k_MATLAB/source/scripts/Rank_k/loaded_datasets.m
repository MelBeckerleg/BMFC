%BMC 
%Simple script to load dataset and test on different algorithms
%Author: Mel
%Date: 03/02/2020


%% Add paths
close('all');clear;clc
%impath='/home/beckerleg/' ;%
%homepath='/home/beckerleg/Ethera/' ;%
impath='/home/user/Documents/Mel/Ethera/FirstYear/Images/' ;%
homepath='/home/user/Documents/Mel/Ethera/' ;%
addpath([homepath 'BMFC/Rank_k_MATLAB/source/utils'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/ADMM/utils/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/ADMM/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/error/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/models/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/nmfv1_4'])


%% Test Parameters
test_params.homepath=homepath;
test_params.impath=impath;
test_params.num_trials=100;





%% Problem Parameters
%dataset_name='reddit';
trial_idx=1;dataset_name='chembl';
m=7500;n=7500;
%m=5000;n=1000;
k=max(floor(sqrt(m)),1);
prob_params.m=m; prob_params.n=n; prob_params.random_subset_trial_idx=trial_idx;
prob_params.rho=0.7;
prob_params.generate=@load_dataset;prob_params.dataset_name=dataset_name;
prob_params.homepath=homepath;
prob_params.k=k;

%% Do Stuff
%Generate problem
num_methods=8;
err_train=zeros([test_params.num_trials,num_methods]);err_test=zeros([test_params.num_trials,num_methods]);err_train_binary=zeros([test_params.num_trials,num_methods]);err_test_binary=zeros([test_params.num_trials,num_methods]);



for trial_idx=1:test_params.num_trials
    prob_params.random_subset_trial_idx=rng(trial_idx);
    prob_params.m=m;prob_params.n=n;
    [db,prob_params]=prob_params.generate(prob_params); 
    %db=db';[prob_params.m,prob_params.n]=size(db);
    [y_train, mask_train, y_val, mask_val, y_test, mask_test] = split_observed(db, [prob_params.rho, (1-prob_params.rho)/2, (1-prob_params.rho)/2],0);
    %if the matrix is sparse, use W_omega=zeros(size(db));
    W_omega=sparse(zeros(size(db)));W_omega(mask_train)=db(mask_train);
    if sum(sum(W_omega))>0
        W_omega=-W_omega;
        db(mask_train)=-db(mask_train);
        db(mask_test)=-db(mask_test);
        db(mask_val)=-db(mask_val);
    end
    method_idx=1;
for solver_id={'BNMF','TBMC','SPLIT-AVG','SPLIT-DIVIDE','TBMC-IU','SPLIT-AVG-IU','SPLIT-DIVIDE-IU','Spectral'}%TBMC-side
    solve_params=configure_solver(solver_id{1});solve_params.k=prob_params.k;

    [X,Y]=solve_params.rank_k_solver(W_omega,solve_params);
    
    
    
    %Evaluate the error
    err_train(trial_idx,method_idx)=binary_rank_error('l2_scaled',mask_train,X,Y,db,[],[],'shift',1);
    err_test(trial_idx,method_idx)=binary_rank_error('l2_scaled',mask_test,X,Y,db,[],[],'shift',1);
    err_train_binary(trial_idx,method_idx)=binary_rank_error('binary',mask_train,X,Y,db,[],[],'shift',1);
    err_test_binary(trial_idx,method_idx)=binary_rank_error('binary',mask_test,X,Y,db,[],[],'shift',1,'bin',1);
    
    sprintf('Method: %s \n Training error from train is %s \n Test error is %s \n Binary Test error is %s',string(solver_id),string(err_train(trial_idx,method_idx)),string(err_test(trial_idx,method_idx)),string(err_test_binary(trial_idx,method_idx)))

    bench_err_train0=binary_rank_error('binary',mask_train,zeros([prob_params.m,1]),zeros([prob_params.n,1]),db,[],[],'shift',1);
    bench_err_train1=binary_rank_error('binary',mask_train,ones([prob_params.m,1]),ones([prob_params.n,1]),db,[],[],'shift',1);
    
    sprintf('\nBenchmark Training error from train is %s or %s',string(bench_err_train1),string(bench_err_train0))
    
    method_idx=method_idx+1;
end

end

if 0
    load(sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/rank_k_%s_m%s_n%s_inv',dataset_name,string(m),string(n)),'err_test','err_train','err_test_binary')
    boxplot(err_test_binary)
    xticklabels([{'BNMF','TBMC','SPLIT-AVG','SPLIT-DIVIDE','TBMC-IU','SPLIT-AVG-IU','SPLIT-DIVIDE-IU','Spectral'}]);xtickangle(15);
    ylabel('Test Error')
    saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/rank_k_%s_m%s_n%s_inv',dataset_name,string(m),string(n)),'epsc')
    save(sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/rank_k_%s_m%s_n%s_inv',dataset_name,string(m),string(n)),'err_test','err_train','err_test_binary')
end


function [solve_params] = configure_solver(solver_id)
solve_params.label=solver_id;
%TMBC,
if strcmp(solver_id, 'TBMC')
    %these are different
    solve_params.rank_k_solver=@splitting;
    solve_params.it=0;
    solve_params.rank_1_method='lp';
    solve_params.pf=0;



end

if strcmp(solver_id, 'TBMC-side')
    %these are different
    solve_params.rank_k_solver=@splitting;
    solve_params.it=0;
    solve_params.rank_1_method='lp';
    solve_params.pf=1;



end


% TBMC with iterative updating
if strcmp(solver_id, 'TBMC-IU')
    %these are different
    solve_params.rank_k_solver=@splitting;
    solve_params.it=1;
    solve_params.rank_1_method='lp';
    solve_params.pf=0;


end

% TBMC with iterative updating
if strcmp(solver_id, 'TBMC-IU-side')
    %these are different
    solve_params.rank_k_solver=@splitting;
    solve_params.it=1;
    solve_params.rank_1_method='lp';
    solve_params.pf=1;


end


if strcmp(solver_id, 'SPLIT-AVG')
%these are different
solve_params.rank_k_solver=@splitting;
solve_params.it=0;
solve_params.rank_1_method='avg';
solve_params.pf=0;



end

if strcmp(solver_id, 'SPLIT-AVG-IU')
    %these are different
    solve_params.rank_k_solver=@splitting;
    solve_params.it=1;
    solve_params.rank_1_method='avg';
    solve_params.pf=0;


end


if strcmp(solver_id, 'SPLIT-DIVIDE')
%these are different
solve_params.rank_k_solver=@splitting;
solve_params.it=0;
solve_params.rank_1_method='partition';
solve_params.pf=0;



end

if strcmp(solver_id, 'SPLIT-DIVIDE-IU')
    %these are different
    solve_params.rank_k_solver=@splitting;
    solve_params.it=1;
    solve_params.rank_1_method='partition';
    solve_params.pf=0;

end


%BNMF

if strcmp(solver_id, 'BNMF')
    %these are different
    solve_params.rank_k_solver=@BNMF;
    %solve_params.it=1;
    %solve_params.rank_1_method='partition';
    


end

if strcmp(solver_id,'Neighbours')
    solve_params.rank_k_solver=@neighbour;
end


if strcmp(solver_id,'Nuclear')
    solve_params.rank_k_solver=@convex;
    recovered=0;
    binarised=1;
end

if strcmp(solver_id,'NuclearR')
    solve_params.rank_k_solver=@convex;
    recovered=1;
end

if strcmp(solver_id,'Spectral')
    solve_params.rank_k_solver=@spectral;
end

end
