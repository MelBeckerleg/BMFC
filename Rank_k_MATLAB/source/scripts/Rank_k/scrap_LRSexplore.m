m=500;n=m;k=0.15*m;tau=0.3;rho=0.1;synthetic_set_up='Brickwork';exact_rho=rho;

rng(1);
prob=rank_k_problem();
prob.m=m;prob.n=n;prob.epsilon=epsilon;prob.rho=rho;prob.tau=tau;prob.k=k;prob.k_solve=1;
prob.set_up=synthetic_set_up; % generate_row_clusters
[prob.X_true,prob.Y_true,prob.A_true,prob.W_omega]= generate_row_clusters(prob);
prob.mask=1.*(abs(prob.W_omega)>0);

opts.epsilon= epsilon;
opts.k = k
%
[Xb,Yb] = BNMF(prob.W_omega,opts);
CompletedMat = X*Y';
% 
% function [X,Y] = BNMF(W_omega,opts)
%     addpath(genpath('/home/user/Documents/Mel/Ethera/Rank_k/Rank_k_MATLAB'));
%     k=opts.k;
%     A=(W_omega+1)/2;
%     A(A==1/2)=NaN;
%     [W,H]=wnmfrule(A,k);H=H';
%     [W,H]=binary_rescale(W,H);
%     X=(W>0.5)*1.;
%     Y=(H>0.5)*1.;
% 
% 
% end
% 

if 0
%distance:
lambda = 2/(m);%/(sqrt(sqrt(exact_rho*m)));%in the paper they take 1/sqrt(rho*n) but this is too small
mu_weight_one = 1;
mu_weight_two = 1;
[CompletedMat, E_hat, iter] = inexact_alm_rpca(prob.W_omega*prob.mask, lambda,mu_weight_one,mu_weight_two,1e-9,100);


%CGIHT
    oopts.verbosity = 0;
    oopts.reltol = 1e-3; 
    oopts.maxiter = 100; 
    oopts.maxit = 100;
    oopts.rate_limit = 1-oopts.reltol; 
    oopts.relres = oopts.reltol*norm(prob.W_omega); 
    oopts.rel_res_tol = 1e-3; 
    oopts.rel_res_change_tol = 1e-3;
    oopts.itres = zeros(oopts.maxiter,1);
    [start.U,sigma,start.V] = svd(prob.W_omega);start.sigma=diag(sigma);
    start.L = start.U;start.R=start.V;
    [CompletedMat1, ~] = NIHT_Matrix(m,n,k,find(prob.W_omega),prob.W_omega(find(prob.W_omega)),start,oopts);
    [CompletedMat2, ~] = ASD(m,n,k,find(prob.W_omega),prob.W_omega(find(prob.W_omega)),start,opts);

    


    


max(CompletedMat)

end
recovered = 1;

sprintf('The error is \n')
val = nnz(CompletedMat1 - prob.A_true)/(n*m)
val = nnz(CompletedMat2 - prob.A_true)/(n*m)

sprintf('The binarised error is \n')
val = nnz((CompletedMat>0.0)*1. - prob.A_true)/(n*m)
        

recover_rows_function = @nearest_group;
X = recover_rows_function((CompletedMat1>0)*1.,k,opts);
Y = majority_vote(prob.A_true,X);
sprintf('The final error is \n')
val = nnz((X*Y'>0)*1. - prob.A_true)/500^2

%Trying to work out why this *s up at 0.15 (and 0.3?)
for ii = [1:k]
ordered = find(prob.X_true(:,ii));
idx(val:val+length(ordered)-1) = ordered;
val = val+length(ordered);
end