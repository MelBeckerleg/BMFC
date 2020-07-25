%Evaluate the cluster strength of a recovered matrix X_approx and a ground truth X
%strength=cluster_strength(X,Xapprox)
%X is the ground truth and Xapprox is the approximation of X.
%method: since different orderings of columns are combinatorial, the method scores based on the best and the worst overlap: X'*Xapprox gives the overlap between columns of X and Xapprox. Ideal clusters will have either very large or very small values. The metric penalises 'middle' values that are 'spread out'. More specifically, if the X_approx gives good recovery of clusters, then each column will have a large overlap with one of the columns of X, and since the clusters of X are assumed to be distinct, the overlap with other columns will be smaller.  Hence, the strength factor is he minimum of
%a) the rating for the least representative column min(sum(abs(max(X'*Xapprox))) of the positive entries and 
%b) the rating for the least representative column min(sum(abs(max((1-X)'*(1-Xapprox)))) of the non-positive entries 
%both of these are normalised by the size of X and the proportion of positive/non-positive values respectively.
%Author: Mel Beckerleg. Date: 17/12/2019

function [val] = cluster_strength(X,Xapprox) 
[m,r]=size(X);
base_val1=max(max(X'*X));
base_val2=max(max((1-X)'*(1-X)));
val1=1/(base_val1*r)*min(sum(abs(max(X'*Xapprox))),sum(abs(max((X'*Xapprox))')));
val2=1/(base_val2*r)*min(sum(abs(max((1-X)'*(1-Xapprox)))),sum(abs(max(((1-X)'*(1-Xapprox)))')));
val=min(val1,val2);
end
