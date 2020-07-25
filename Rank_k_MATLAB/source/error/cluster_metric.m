function [funs]=cluster_metric()
funs.cluster_error=@cluster_error;
funs.cluster_error_small_r=@cluster_error_small_r;
end

function [val] = cluster_error_small_r(X,Xapprox,r) 
val=0;
for rval=1:r
val=val+mean(pdist(X(Xapprox(:,rval)>0,:)));
end
val=1/r*val;
end


function [val] = cluster_error(X,Xapprox) 
[m,r]=size(X);
base_val1=max(max(X'*X));
base_val2=max(max((1-X)'*(1-X)));
val1=1/(base_val1*r)*min(sum(abs(max(X'*Xapprox))),sum(abs(max((X'*Xapprox))')));
val2=1/(base_val2*r)*min(sum(abs(max((1-X)'*(1-Xapprox)))),sum(abs(max(((1-X)'*(1-Xapprox)))')));
val=min(val1,val2);
end