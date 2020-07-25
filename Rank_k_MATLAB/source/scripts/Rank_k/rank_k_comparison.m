%Script to run solve the rank-k problem.
clear; 
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/'))
addpath(genpath('/home/user/mosek/'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/solvers/MatrixCompletion'))
tag='all_july_24'

loaded=0;generate=1;
num_trials=30;m_vals=[50,100,150,200,300,400,500];%,400,500];%[50,100,150,200,250,300];
epsilon=0.03;tau=0.3;rho=0.7;


solvers={'splitting','Neighbour','Recover','BNMF'};
splitting_methods={'lp','avg','partition','nmf'};
recovery_methods={'convex','1bitMC','RPCA'};%CGIHT is unconfigured
num_methods=13;

%just testing the newbies
%num_methods=2;solvers={'Recover'};recovery_methods={'RPCA','1bitMC'};


names={};


err_train=zeros(num_trials,num_methods);
err_test=zeros(num_trials,num_methods);
row_cluster_strength=zeros(num_trials,num_methods);
col_cluster_strength=zeros(num_trials,num_methods);
err_hamming=zeros(num_trials,num_methods);
comp_time=zeros(num_trials,num_methods);
for m = m_vals
    n=m;k=floor(log(m));ksolve=k;opts=rank_k_solve_opts;opts.k=ksolve;opts.pf=0;
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
        if loaded
            prob=rank_k_problem(); 
            prob.set_up='Movie_Lens';
            [prob.A_true,prob.W_omega]= ML_load(prob);
        end
    end
    method_idx=1;
    for my_solver = solvers
        if strcmp(my_solver,'splitting')
            for method=splitting_methods
                opts.it=0;opts.rank_1_method=method;
                tic;
                [X,Y]=splitting(prob.W_omega,opts); %cycle through ip,lp,avg,partition
                %toc;
                comp_time(trial_idx,method_idx)=toc;
                [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
                names{method_idx}=[my_solver,'_',method,'_','_no_it'];method_idx=method_idx+1;
                opts.it=1;
                tic;
                [X,Y]=splitting(prob.W_omega,opts); 
                %toc;
                toc;comp_time(trial_idx,method_idx)=toc;
                [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
                names{method_idx}=[my_solver{1},'_',method,'_','_it'];method_idx=method_idx+1;
            end
        else 
            if strcmp(my_solver,'Recover')
                for method = recovery_methods
                    opts.recovery_solver=method;
                    func = str2func(my_solver{1});
                    tic;
                    [X,Y] = func(prob.W_omega,opts);
                    %toc;
                    comp_time(trial_idx,method_idx)=toc;%
                    [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
                    names{method_idx}=[my_solver{1},'_',method];method_idx=method_idx+1;

                end
                
            
        else
            func=str2func(my_solver{1});
            tic;
            [X,Y]=func(prob.W_omega,opts);  
            comp_time(trial_idx,method_idx)=toc;%
         	[err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming]=error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming);
            names{method_idx}={my_solver{1}};method_idx=method_idx+1;
            end
        end
    end
fn=sprintf('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_%s_k_%s_rho_%s_tau_%s_eps_%s_%s.mat',string(m),string(k),string(100*rho),string(100*tau),string(100*epsilon),tag);
save(fn,'names','err_train','err_test','row_cluster_strength','col_cluster_strength','comp_time');

end
%clear prob; 
end