%Evaluate the error of the clusterings and the recovery
%Author: Mel Beckerleg
%Date 28/10/2019

function [val] = cluster_error(A,mask,X,Y,X_approx,Y_approx)
[m,n]=size(A);rho=nnz(mask)/m/n;
val.error=(norm(A-X_approx*Y_approx','fro')/sqrt(m*n))^2;
val.xerror=cluster_strength(X,X_approx);
val.yerror=cluster_strength(Y,Y_approx);
pred_error=A-X_approx*Y_approx';
if sum(sum(mask))
    val.train_error=(norm(pred_error((mask>0)),'fro')/(sqrt((1-rho)*m*n)))^2;
end
if sum(sum(1-mask))
    preds=A-X_approx*Y_approx';
    val.test_error=(norm(pred_error(((1-mask)>0)),'fro')/(sqrt(rho*m*n)))^2;
end
val.recovered_rank='FixThis';

end


function [val] = cluster_strength(X,Xapprox) 
[m,r]=size(X);
base_val1=max(max(X'*X));
base_val2=max(max((1-X)'*(1-X)));
val1=1/(base_val1*r)*min(sum(abs(max(X'*Xapprox))),sum(abs(max((X'*Xapprox))')));
val2=1/(base_val2*r)*min(sum(abs(max((1-X)'*(1-Xapprox)))),sum(abs(max(((1-X)'*(1-Xapprox)))')));
val=min(val1,val2);
end