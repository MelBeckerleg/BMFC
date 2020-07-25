function [val] = cluster_strength(X,Xapprox) 
[m,r]=size(X);
base_val1=max(max(X'*X));
base_val2=max(max((1-X)'*(1-X)));
val1=1/(base_val1*r)*min(sum(abs(max(X'*Xapprox))),sum(abs(max((X'*Xapprox))')));
val2=1/(base_val2*r)*min(sum(abs(max((1-X)'*(1-Xapprox)))),sum(abs(max(((1-X)'*(1-Xapprox)))')));
val=min(val1,val2);
end
