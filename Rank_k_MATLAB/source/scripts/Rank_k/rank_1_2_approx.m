%Script to run solve the rank-1 problem.
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



num_trials=100;
m=10;n=10;epsilon=0.03;tau=0.7;rho=0.7;
methods={'lp','avg','partition','nmf'};
%errors=['l2'];
error=zeros(num_trials,2*length(methods));
cont=1;
for trial_idx=1:num_trials
    if cont
        rng(trial_idx)
        %generate
        %prob=rank_k_problem();
        prob.m=m;prob.n=n;prob.epsilon=epsilon;prob.rho=rho;prob.tau=tau;prob.k=3;prob.k_solve=1;
        prob.set_up='generate_block_diagonal';%'generate_row_clusters';%random';%'Brickwork'; % generate_row_clusters
        %[prob.A_true,prob.W_omega]= generate_random(prob);        
        [prob.X_true,prob.Y_true,prob.W_omega,prob.W_omega] = generate_row_clusters(prob); %changed to W_omega (when there is noise) for 2approx
        %prob.a=0.5;[prob.X_true,prob.Y_true,prob.A_true,prob.W_omega] = generate_block_diagonal_clusters_geometric(prob);
        prob.mask=1.*(abs(prob.W_omega)>0);


        %solve with rank_1
        [X,Y]=rank_1_solve(prob.W_omega,'ip');
        base_err=binary_rank_error('l2',prob.mask,X,Y,prob.W_omega>0)^2;%changed to W_omega (when there is noise) for 2approx
        method_idx=0;

        for method=methods
            method_idx=method_idx+1;
            [X,Y]=rank_1_solve(prob.W_omega,method); %cycle through ip,lp,avg,partition
            err=binary_rank_error('l2',prob.mask,X,Y,prob.W_omega>0)^2; %changed to W_omega (when there is noise) for 2approx
            [X,Y]=iterative_update(X,Y,prob.W_omega);
            err_it=binary_rank_error('l2',prob.mask,X,Y,prob.W_omega>0)^2; 

            if base_err
                if err_it/base_err<0.999
                    disp(method)
                    sprintf('%d',err_it)
                    sprintf('%d',base_err)
                    cont=0;
                    break
                end
                
                if err>2*base_err && strcmp(method,'lp')
                    disp(method)
                    sprintf('%d',err)
                    sprintf('%d',base_err)
                    cont=0;
                    break
                end
                
                error(trial_idx,method_idx)=err/base_err;
                error(trial_idx,method_idx+length(methods))=err_it/base_err;
            else
                if err
                    error(trial_idx,method_idx)=m*n+1;
                end
                if err_it
                    error(trial_idx,method_idx+length(methods))=m*n+1;
                end
            end

        end
    else
        break
    end
end


%mean,median,proportion that don't find a zero error
%solution when there is one, proportion that vioalte the 2-approximation,
%time taken
chart=zeros(4,8);
for method_idx=1:8
    non_zero_idx=find(error(:,1)~=0);
    chart(1,method_idx)=mean(error(non_zero_idx,method_idx));
    
   chart(2,method_idx)=median(error(non_zero_idx,method_idx));
    zero_idx=find(error(:,1)==0);
   chart(3,method_idx)=sum(error(zero_idx,method_idx)>0)/length(zero_idx);
   chart(4,method_idx)=sum(error(:,method_idx)>2.0001)/length(error(:,method_idx));
end


chart

