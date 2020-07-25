%Generate a database according to a row cluster (brickwork) model. 
function [X_true,Y_true,A_true,W_omega] = generate_row_clusters(opts)
%Set_up and ground truth
m=opts.m;n=opts.n;k=opts.k;tau=opts.tau;


X_true=zeros(m,k);%Y_true=zeros(n,k);
x_vals=randperm(m);%y_vals=randperm(n);
for kk=1:k
    X_true(x_vals(1+(kk-1)*floor(m/k):kk*floor(m/k)),kk)=1;
    %Y_true(y_vals(1+(kk-1)*floor(n/k):kk*floor(n/k)),kk)=1;
end
Y_true=rand(n,k)>tau;
A_true=X_true*Y_true';

%apply noise and subsampling filter
epsilon=opts.epsilon;rho=opts.rho;
A = A_true + (1-2*A_true).*(rand(m,n)<epsilon); 
W_omega = 2*A-1;
W_omega(rand(m,n)<rho)=0;
end







