function [X,Y]=Recover(W_omega,opts)
k=opts.k;recovered=opts.recovered;binarised=opts.binarised;
recovery_solver=opts.recovery_solver;%'convex';
[m,n]=size(W_omega);

%set default methods for cluster recovery
binarised=1; binarise_threshold = 0.0;%results are for clustering
recovered = 0; %best to use the observed database over the recovered one...
                %with the exception of LRS? 
recover_rows_function = @nearest_group;
recover_footprints_function = @majority_vote;


B=abs(W_omega); %B=mask
exact_rho=sum(sum(B))/m/n;
%A=(W_omega+1)/2; %currently just using W_omega rather than A, some
%indication this is a better shout



if strcmp(recovery_solver,'convex')
    [CompletedMat, ~] = MatrixCompletion(W_omega.*B, B,100, 'nuclear', 10, 1e-8, 0);

else
    if strcmp(recovery_solver,'CGIHT')
    reltol = 1e-3; 
    maxiter = 100; 
    rate_limit = 1-reltol; 
    relres = reltol*norm(W_omega); 
    itres = zeros(maxiter,1);
    print('not configured for CGIHT')
    %[CompletedMat, ~] = Modified_CGIHT_Matrix(m,n,k,B,W_omega,start,opts);
    
    else
        if strcmp(recovery_solver, '1bitMC')
            %opts.epsilon = 0.2; %Ideally you want to work out if you're passing 
                                %this as a parameter or not
            
            CompletedMat = OneBitMC(m,n,k,[],W_omega,opts);
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
            opts.thresh=sqrt(n^(3/4-1/6+1/24)); %for smaller m, use this as a proxy
            %opts.thresh=sqrt(2*sqrt(opts.thresh));
            recover_rows_function = @threshold_group;
            binarise_threshold = 0.0;
            
            
%             sprintf('The error is \n')
%             prob.A_true=varagin;
%             val = nnz(CompletedMat - prob.A_true)/(n*m)
% 
% 
%             sprintf('The binarised error is \n')
%             val = nnz((CompletedMat>0)*1. - prob.A_true)/(n*m)
        else
            if strcmp(recovery_solver, '1bitMCnearest')
                %opts.epsilon = 0.2; %Ideally you want to work out if you're passing 
                                    %this as a parameter or not

                CompletedMat = OneBitMC(m,n,k,[],W_omega,opts);  
                binarise_threshold = 0.0;
            else
                if strcmp(recovery_solver,'RPCA')
                    %Uses  RPCA inexact ALM method by Chen 2009
                    % [A_hat E_hat iter] = inexact_alm_rpca(D, lambda, tol, maxIter)
                    % we choose lambda as in wright to be 1/sqrt(m);
                    lambda = 1/sqrt(sqrt(exact_rho*m));%in the paper they take 1/sqrt(rho*n) but this is too small
                    [CompletedMat, ~, ~] = inexact_alm_rpca(W_omega.*B, lambda,1e-8,100);
                    recovered = 1;
                else
                    if strcmp(recovery_solver,'NIHT')
                            oopts.verbosity = 0;
                            oopts.reltol = 1e-3; 
                            oopts.maxiter = 100; 
                            oopts.maxit = 100;
                            oopts.rate_limit = 1-oopts.reltol; 
                            oopts.relres = oopts.reltol*norm(W_omega); 
                            oopts.rel_res_tol = 1e-3; 
                            oopts.rel_res_change_tol = 1e-3;
                            oopts.itres = zeros(oopts.maxiter,1);
                            [start.U,sigma,start.V] = svd(W_omega);start.sigma=diag(sigma);
                            start.L = start.U;start.R=start.V;
                            [CompletedMat, ~] = NIHT_Matrix(m,n,k,find(W_omega),W_omega(find(W_omega)),start,oopts);
                            recovered = 1; binarised = 1;binarise_threshold = 0 ;

                    else
                        if strcmp(recovery_solver,'ASD')
                            oopts.verbosity = 0;
                            oopts.reltol = 1e-3; 
                            oopts.maxiter = 100; 
                            oopts.maxit = 100;
                            oopts.rate_limit = 1-oopts.reltol; 
                            oopts.relres = oopts.reltol*norm(W_omega); 
                            oopts.rel_res_tol = 1e-3; 
                            oopts.rel_res_change_tol = 1e-3;
                            oopts.itres = zeros(oopts.maxiter,1);
                            [start.U,sigma,start.V] = svd(W_omega);start.sigma=diag(sigma);
                            start.L = start.U;start.R=start.V;
                            [CompletedMat, ~] = ScaledASD(m,n,k,find(W_omega),W_omega(find(W_omega)),start,oopts);
                            recovered = 1;
                            binarised = 1; binarise_threshold = 0;

                        else
                            disp('method not recognised, returning observed matrix')
                            CompletedMat = W_omega;
                        end
                    end
                    end
            end
        end
    end
end    


if binarised
    CompletedMat=(CompletedMat>binarise_threshold)*1;
end

if recovered
    mat=CompletedMat;
else
    mat=W_omega;    
end


X = recover_rows_function(CompletedMat,k,opts);
Y = recover_footprints_function(mat,X);
    
end
    
