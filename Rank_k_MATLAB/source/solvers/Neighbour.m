function [X,Y]= Neighbour(W_omega,opts)
k=opts.k;
X=nearest_group(W_omega,k);
Y = majority_vote(W_omega,X); 
end
