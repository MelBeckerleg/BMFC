function[W,H]= binary_rescale(W,H)
%check that there are no empty columns (to facilitate inv!)
    keep_idx=find(sum(abs(W))>0);
    W=W(:,keep_idx);H=H(:,keep_idx);
    keep_idx=find(sum(abs(H))>0);
    W=W(:,keep_idx);H=H(:,keep_idx);

%works for m,n>1!
    [m,k]=size(W);[n,k]=size(H);
    Dw=diag(max(W));
    Dwinv=diag(1./max(W));
    Dh=diag(max(H));
    Dhinv=diag(1./max(H));
    W=W*(sqrt(Dwinv)*sqrt(Dh)')';
    H=H*(sqrt(Dw)*sqrt(Dhinv)')';
end
