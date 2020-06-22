%Script to run solve the rank-k problem.
clear; 
addpath(genpath('/home/user/Documents/Mel/Ethera/Rank_k/Rank_k_MATLAB/'))




load=0;generate=1;
num_trials=100;
m=50;n=50;epsilon=0.0;tau=0.5;rho=0.7;k=3;ksolve=k;
methods={'lp','avg','partition','nmf'};
solvers={'splitting','Neighbour','Nuclear','BNMF'};
names={};
opts=rank_k_solve_opts;opts.k=ksolve;

err_train=zeros(num_trials,length(methods)+length(solvers)-1);
err_test=zeros(num_trials,length(methods)+length(solvers)-1);
row_cluster_strength=zeros(num_trials,2*length(methods)+length(solvers)-1);
col_cluster_strength=zeros(num_trials,length(methods)+length(solvers)-1);
err_hamming=zeros(num_trials,length(methods)+length(solvers)-1);


for trial_idx=1:num_trials
    rng(trial_idx)
    %generate
    if generate
        prob=rank_k_problem();
        prob.m=m;prob.n=n;prob.epsilon=epsilon;prob.rho=rho;prob.tau=tau;prob.k=k;prob.k_solve=1;
        prob.set_up='Brickwork'; % generate_row_clusters
        [prob.X_true,prob.Y_true,prob.A_true,prob.W_omega]= generate_row_clusters(prob);
        prob.mask=1.*(abs(prob.W_omega)>0);
        if ~sum(abs(prob.W_omega))
            disp('cant tile empty database')
        end
    else
        if load
            prob=rank_k_problem();
            prob.set_up='Movie_Lens';
            [prob.A_true,prob.W_omega]= ML_load(prob);
        end
    end
    method_idx=1;
    for my_solver = solvers
        if strcmp(my_solver,'splitting')
            for method=methods
                opts.it=0;opts.rank_1_method=method;
                [X,Y]=splitting(prob.W_omega,opts); %cycle through ip,lp,avg,partition
                [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
                names{method_idx}=[my_solver,'_',method,'_','_it'];method_idx=method_idx+1;
                opts.it=1;
                [X,Y]=splitting(prob.W_omega,opts); 
                [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
                names{method_idx}=[my_solver{1},'_',method,'_','_it'];method_idx=method_idx+1;
            end
        else
            func=str2func(my_solver{1});
            [X,Y]=func(prob.W_omega,opts);  
            [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
            names{method_idx}={names, my_solver};method_idx=method_idx+1;
        end
    end

end
