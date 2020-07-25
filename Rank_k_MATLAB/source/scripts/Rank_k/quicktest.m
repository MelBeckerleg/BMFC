%% Quick test of LRS
% Easily test one algorithm

%% Set-up the model
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/'))
addpath(genpath('/home/user/mosek/'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/solvers/MatrixCompletion'))

% params
epsilon=0.03;tau=0.3;rho=0.7;m=300;k=floor(log(m));n=m;
% generate
prob=rank_k_problem();
prob.m=m;prob.n=n;prob.epsilon=epsilon;prob.rho=rho;prob.tau=tau;prob.k=k;prob.k_solve=1;
prob.set_up='Brickwork'; % generate_row_clusters
[prob.X_true,prob.Y_true,prob.A_true,prob.W_omega]= generate_row_clusters(prob);
prob.mask=1.*(abs(prob.W_omega)>0);
B=prob.mask;A=(prob.W_omega+1)/2;

%% Solve with LRS
% Description of second code block
%[X,Y]=Recover(W_omega,opts)

%LRS: doesn't seem to work too well... has about a 70% error

lambda = 1/sqrt(sqrt(rho*n));%in the paper they take 1/sqrt(rho*n) but this is too small
[CompletedMat, S_0, ~] = inexact_alm_rpca(prob.W_omega.*B, lambda,1e-8,1000);
%check whether we expect to recover
%i.e. is k<(epsilon+rho)*m*n/m/n*m*log(n)^2/mu


%% Solve with 1bitMC
opts.epsilon = 0.2; %Ideally you want to work out if you're passing 
                    %this as a parameter or not
exact_rho=nnz(W_omega)/m/n;
CompletedMat = OneBitMC(m,n,k,B,prob.W_omega,opts);
%set the threshold according to the definition in chapter 5,
c_rho = n^(1/24)*exact_rho;
c_k = n^(1/12)*k;
opts.thresh = sqrt(32*sqrt(2)*sqrt(c_k/c_rho)*(1+1/(4*(1/2-opts.epsilon)^2)) ...
                           *(1/2-opts.epsilon)*(1+(1/2-opts.epsilon)^2/(1/4-(1/2-opts.epsilon)^2)) ...
                           *8*sqrt(2)*(1+sqrt(6))/exp(1) ...
                           *sqrt(1+(m+n)*log(m)/(exact_rho*m*n)) ...
                           *n^(3/4-1/6+1/24)); %take sqrt(T_epsilon) since pdist calculates 
opts.thresh=sqrt(n^(3/4-1/6+1/24));

%% Recover clusterings
opts.k=k; %just a hack

recover_rows_function = @nearest_group;
recover_rows_function = @threshold_group; %uncomment fo 1-bit
recover_footprints_function = @majority_vote;

X = recover_rows_function(CompletedMat,k,opts);
%Y = recover_footprints_function(prob.W_omega,X); %default
Y = recover_footprints_function(CompletedMat,X); %uncomment if looking at LRS

    