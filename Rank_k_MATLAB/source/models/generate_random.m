function [A_true,W_omega] = generate_random(opts)
m=opts.m;n=opts.n;k=opts.k;tau=opts.tau;rho=opts.rho;

A_true=(rand(m,n)>tau)*1.;
W_omega = 2*A_true-1;
W_omega(rand(m,n)<rho)=0;
end