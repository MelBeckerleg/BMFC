function [X,Y] = BNMF(W_omega,opts)
    addpath(genpath('/home/user/Documents/Mel/Ethera/Rank_k/Rank_k_MATLAB'));
    k=opts.k;
    A=(W_omega+1)/2;
    A(A==1/2)=NaN;
    [W,H]=wnmfrule(A,k);H=H';
    [W,H]=binary_rescale(W,H);
    X=(W>0.5)*1.;
    Y=(H>0.5)*1.;


end