%MC_DEMO_SOLVE_ADMM     Solve matrix completion on graphs with ADMM
%
%
%
% see also: (/source/features/utils) MC_solve_ADMM, MC_demo_grid_search, split_observed,
%           sample_sparse, sample_sparse_t, sample_sparse_AtA, 
%            (/utils) vec
%
%
%code author: Vassilis Kalofolias
%date: Nov 2014
%code edits: Mel Beckerleg
%date: August 2019
%%

close('all');clear;clc
%impath='/home/beckerleg/' ;%
%homepath='/home/beckerleg/Ethera/' ;%
impath='/home/user/Documents/Mel/Ethera/FirstYear/Images/' ;%
homepath='/home/user/Documents/Mel/Ethera/' ;%
addpath([homepath 'BMFC/Rank_k_MATLAB/source/utils'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/ADMM/utils/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/ADMM/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/error/'])
%% test parameters 
test_params.homepath=homepath;
test_params.impath=impath;
test_params.row_features={'ch_fingerprints','none'};
test_params.col_features={'ch_protein_cathids','ch_protein_funfams','ch_protein_sigs','none'};
test_params.num_trials=5;

%% Initialise error vectors

rmse_test_vals=zeros(2,length(test_params.row_features),length(test_params.col_features),test_params.num_trials);

%% load Chembl data
load([homepath 'MP2/Data/read_ch/mats/ch_dense_set'],'supported_ch','ch_fingerprints','ch_protein_funfams','ch_protein_sigs'); 
[M,N]=size(supported_ch);none_row=sparse(M,1);none_col=sparse(N,2);load([homepath 'MP2/Data/read_ch/mats/ch_protein_cathids.txt']);

%% problem parameters 

for update_gamma_val=[0,1]
for trial_idx=1:test_params.num_trials    
    for row_feature_idx=1:length(test_params.row_features)
        for col_feature_idx=1:length(test_params.col_features)
        
            prob_params.m=500; prob_params.n=500; prob_params.random_subset_trial_idx=trial_idx;
            prob_params.row_features=test_params.row_features{row_feature_idx};
            prob_params.col_features=test_params.col_features{col_feature_idx};
            prob_params.rho=0.7;
            %% Generate the graphs and the subset of rows and columns
            [Xn,prob_params.Lr,prob_params.Lc]=generate_sub(eval(prob_params.row_features),eval(prob_params.col_features),supported_ch,prob_params);
            prob_params.size_X=size(Xn);

            %% Subsample
            [y_train, mask_train, y_val, mask_val, y_test, mask_test] = split_observed(Xn, [prob_params.rho, (1-prob_params.rho)/2, (1-prob_params.rho)/2],0);
            y_train=(y_train+1)/2;y_test=(y_test+1)/2;y_val=(y_val+1)/2;
            Xn((Xn==-1))=0;


            %normalize data to zero mean and keep the linear transformation details
            y_lims_init = [min(y_train), max(y_train)];
            if ~isfield(prob_params, 'subtract_mean')
                mean_train=0;
            else
                if prob_params.subtract_mean
                    mean_train = mean(y_train);
                end
            end
            y_train = y_train - mean_train; y_val = y_val - mean_train; %y_test = y_test - mean_train;
            y_lims_scaled = [min(y_train), max(y_train)];

            prob_params.mask_val = mask_val;prob_params.mask_test = mask_test;prob_params.mask_train=mask_train; %added Mel Beckerleg 13/11/2019
            prob_params.A_op = @(x) sample_sparse(x, mask_train);prob_params.At_op = @(x) sample_sparse_t(x, mask_train);prob_params.AtA_op = @(x) sample_sparse_AtA(x, mask_train);

            %% Set problem weighting

            if strcmp(prob_params.row_features,'none')
                prob_params.gamma_r = 0.0;%0.01; %fingerprints is good at about  1/10, generally works well 
            else                
                %c_r=1/(10*rho);prob_params.gamma_r = c_r*rho/2;%0.05;%0.05;%0.01;
                prob_params.gamma_r=0.05;
            end

            if strcmp(prob_params.col_features,'none')    
                prob_params.gamma_c = 0.0;%0.01;
            else            
                %c_c=1/(10*rho);prob_params.gamma_c = c_c*rho/2;%0.05;%0.01;%0.01; 
                prob_params.gamma_c=0.05;
            end
            prob_params.gamma_n=0.1;



            %% SOLVER PARAMETERS 
            solver_params.maxit = 200;
            solver_params.verbose = 3;

            solver_params.tol_abs = 2e-6;
            solver_params.tol_rel = 1e-5;

            % need the scaling used for preprocessing to calculate error correctly
            solver_params.y_lims_init = y_lims_init; solver_params.y_lims_scaled = y_lims_scaled;

            % for small matrices use false!
            solver_params.svds = false;
            
            solver_params.update_gamma=update_gamma_val;

            % MOST IMPORTANT: use verbose = 1 to set rho accordingly (depends on tolerances)
            %solver_params.rho_ADMM = .005000;
            solver_params.rho_ADMM = .2 * geomean([max(1e-3,prob_params.gamma_n), geomean([max(1e-3,norm(y_train)), max(1e-3,prob_params.gamma_r), max(1e-3,prob_params.gamma_c)])]);






            [X_MC_graphs, stat_MC_graphs] = MC_solve_ADMM(y_train, y_val, y_test, prob_params, solver_params);

            %save one example so we know what's going on (roughly)

            fn=[homepath sprintf( 'stat_MC_regularised_cmp%s_prot%s_small_gammn%s_gamr%s_gamc%s_rho%s_m%s_n%s.mat',string(prob_params.col_features),string(prob_params.row_features),...
                string(prob_params.gamma_n),string(prob_params.gamma_r),string(prob_params.gamma_c),string(solver_params.rho_ADMM),string(prob_params.m),string(prob_params.n))];
            save(fn,'X_MC_graphs', 'stat_MC_graphs')

            % calculate the mean
            rmse_test_vals(update_gamma_val+1,row_feature_idx,col_feature_idx,trial_idx)=stat_MC_graphs.rmse_test;

        end
    end
end
end

fn=[homepath sprintf( 'rmse_test_MC_regularised_cmp%s_prot%s_small_gammn%s_gamr%s_gamc%s_rho%s_m%s_n%s.mat',string(prob_params.col_features),string(prob_params.row_features),...
    string(prob_params.gamma_n),string(prob_params.gamma_r),string(prob_params.gamma_c),string(solver_params.rho_ADMM),string(prob_params.m),string(prob_params.n))];
save(fn,'rmse_test_vals')


function [Xn,Lc,Lr] = generate_sub(ft_row,ft_col,value_mat,prob_params)
    [M,N]=size(value_mat);m=prob_params.m;n=prob_params.n;
    if ~isfield(prob_params,'random_subset_trial_idx')
         m_idx=[1:m];n_idx=[1:n];
    else
        rng(prob_params.random_subset_trial_idx)
        m_idx=randi(M,[m,1]);n_idx=randi(N,[n,1]);
        
    end
    
    %overlap between different feature vectors for each row
    Gm=ft_row(m_idx,:)*ft_row(m_idx,:)';
    %ensure that for empty feature vectors, the overlap with themselves is
    %still 1! (really shouldn't have empty feature vectors though!)
    Gm(sum(Gm)==0,sum(Gm)==0)=1;
    Gm=Gm./max(Gm);
    
    %overlap between different feature vectors for each col
    Gu=ft_col(n_idx,:)*ft_col(n_idx,:)';
    Gu=Gu./max(Gu);


    Xn=value_mat(m_idx,n_idx);
    %laplacian= D-A
    Gr_L=diag(sum(Gu))-Gu;Gr_lmax=max(vec(Gu));
    Gc_L=diag(sum(Gm))-Gm;Gc_lmax=max(vec(Gm));
    %normalised Laplacian
    if ~isfield(prob_params,'normalise_laplacian') 
        Lc = (single(full(Gc_L))/Gc_lmax);
        Lr = (single(full(Gr_L))/Gr_lmax);    
    else
        if prob_params.normalise_laplacian
            Lc = (single(full(Gc_L))/Gc_lmax);
            Lr = (single(full(Gr_L))/Gr_lmax);   
        end
    end
end



