%Evaluate the cluster strength of a recovered matrix X_approx and a ground truth X
%strength=cluster_strength(X,Xapprox)
%X is the ground truth and Xapprox is the approximation of X.
%we calculate <X,Xapprox>/||X||/||Xapprox|| for ||.|| the frobenius norm and <.> the inner product induced by the trace norm: 
% <X,Xapprox> = tr(Xapprox'*X)= trace((Xapprox'*X)'(Xapprox'*X)) 
%	      =||(Xapprox'*X)'(Xapprox'*X)  ||_F^2


function [val] = cluster_strength(X,Xapprox) 
[m,r1]=size(Xapprox);
if r1<2
disp('not configured for rank-1, please use hamming distance instead')
val=0;
end

[m,r2]=size(Xapprox);
if r2==0
	val=0;
	else if r2==1
		Xapprox=[Xapprox zeros(m,1)];
		else 
		X=X*1.;Xapprox=Xapprox*1.;
		val1=eval_strength(X,Xapprox);
		val2=eval_strength(1-X,1-Xapprox);

		val=min(val1,val2);
		end
end

end

function out = eval_strength(X,Xapprox)


out=0;
%normalise the columns of Xapprox... 
%...so that the rv correlation makes sense (otherwise we would need to subtract the norm from the 
non_empty_clusters=find(sum(Xapprox));
Xapproxsub=Xapprox(:,non_empty_clusters);
Xapprox(:,non_empty_clusters)=normalize_columns(Xapproxsub);

%and I guess I should for X too
non_empty_clusters=find(sum(X));
Xsub=X(:,non_empty_clusters);
X(:,non_empty_clusters)=normalize_columns(Xsub);

%flip this for consistency with the Robert paper
X=X';Xapprox=Xapprox';


if(sum(sum(abs(Xapprox)))&sum(sum(abs(X))))
	out = trace(X'*X*Xapprox'*Xapprox)/sqrt(trace(X*X')*trace(Xapprox*Xapprox'));
	%val= norm(X'*Xapprox,'fro')^2/norm(X,'fro')/norm(Xapprox,'fro');
end


end

%X'*Xapprox*Xapprox'*X

%method: since different orderings of columns are combinatorial, the method scores based on the best and the worst overlap: X'*Xapprox gives the overlap between columns of X and Xapprox. Ideal clusters will have either very large or very small values. The metric penalises 'middle' values that are 'spread out'. More specifically, if the X_approx gives good recovery of clusters, then each column will have a large overlap with one of the columns of X, and since the clusters of X are assumed to be distinct, the overlap with other columns will be smaller.  Hence, the strength factor is he minimum of
%a) the rating for the least representative column min(sum(abs(max(X'*Xapprox))) of the positive entries and 
%b) the rating for the least representative column min(sum(abs(max((1-X)'*(1-Xapprox)))) of the non-positive entries 
%both of these are normalised by the size of X and the proportion of positive/non-positive values respectively.
%Author: Mel Beckerleg. Date: 17/12/2019
%old code:

%base_val1=max(max(X'*X));
%base_val2=max(max((1-X)'*(1-X)));
%val1=1/(base_val1*r)*min(sum(abs(max(X'*Xapprox))),sum(abs(max((X'*Xapprox))')));
%val2=1/(base_val2*r)*min(sum(abs(max((1-X)'*(1-Xapprox)))),sum(abs(max(((1-X)'*(1-Xapprox)))')));
%val=min(val1,val2);
