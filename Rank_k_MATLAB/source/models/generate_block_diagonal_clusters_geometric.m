%Generate a database according to a block diagonal geometrically decreasing model.
%where tau_{l+1}=a\tau_l 
function [X_true,Y_true,A_true,W_omega] = generate_block_diagonal_clusters_geometric(opts)
%Set_up and ground truth
m=opts.m;n=opts.n;k=opts.k;a=opts.a;
if a==1
    tau=1/k;
else
    tau=(1-a)/(1-a^(k));
end 

X_true=zeros(m,k);Y_true=zeros(n,k);
x_vals=randperm(m);y_vals=randperm(n);
%uncomment if you want to visualise the clusters
%x_vals=1:m;y_vals=1:n;
last=1;

for kk=1:k
    %there are floor(n*tau) positives in this cluster
    new_last=floor(last+floor(n*tau));
    X_true(x_vals(last:min(new_last,n)),kk)=1;
    Y_true(y_vals(last:min(new_last,n)),kk)=1;
    tau=a*tau;
    last=new_last+1;
end

A_true=X_true*Y_true';

%apply noise and subsampling filter
epsilon=opts.epsilon;rho=opts.rho;
A = A_true + (1-2*A_true).*(rand(m,n)<epsilon); 
W_omega = 2*A-1;
W_omega(rand(m,n)<rho)=0;
end
