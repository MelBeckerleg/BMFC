function [X,Y]=Nuclear(W_omega,opts)
k=opts.k;recovered=opts.recovered;binarised=opts.binarised;
recovery_solver='convex';
B=1.*(W_omega==0); %B=mask
A=(W_omega+1)/2;
[m,n]=size(A);
exact_rho=sum(sum(B))/m/n;
X=zeros(m,k);
Y=zeros(n,k);

if strcmp(recovery_solver,'convex')

    [CompletedMat, ~] = MatrixCompletion(A.*B, B,100, 'nuclear', 10, 1e-8, 0);
else
    if strcmp(recovery_solver,'CGIHT')
    reltol = opts.rel_res_tol; 
    maxiter = opts.maxit; 
    rate_limit = 1-opts.rel_res_change_tol; 
    relres = reltol*norm(data); 
    itres = zeros(maxiter,1);
    [CompletedMat, ~] = Modified_CGIHT_Matrix(m,n,r,B,W_omega,start,opts);
    end
end    


if binarised
    CompletedMat=(CompletedMat>0.5)*1;
end
H=squareform(pdist(CompletedMat));
labels=zeros(m,1);
label=1;
%find X

for i=1:m
    if labels(i)==0
        [val,idx]=sort(H(i,:));
        idx=idx(1:ceil(m/k/2));
        if sum(labels(idx))
            existing_label=nonzeros(labels(idx));
            labels(idx)=existing_label(1);
        else
            labels(idx)=label;
            label=label+1;
        end
    end
end
unique_labels=unique(labels);
recovered_rank=min(length(unique_labels),k);
for j=1:recovered_rank
    X(find(labels==unique_labels(j)),j)=1;
end

%find Y

if recovered
    mat=CompletedMat;
else
    mat=A.*B;    
end
for j=1:recovered_rank
    mu=(1/(exact_rho*(length(find(labels==unique_labels(j)))))*sum(mat(find(labels==unique_labels(j)),:)))'>1/2;
    Y(:,j)=mu;
end
    
end
    
