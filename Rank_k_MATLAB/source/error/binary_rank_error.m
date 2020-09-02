function [err]=binary_rank_error(error_type,val_idx,mask,X,Y,varargin)
%structure of varargin: (db,X_true,Y_true,'shift',shift_val)
%you need at least one argument here


if length(varargin )>1 
    if sum(size(varargin{2})) && sum(size(varargin{3}))
        X_true=varargin{2};
        Y_true=varargin{3};
        trues=X_true*Y_true';
    else
        trues=varargin{1};
    end
else 
    trues=varargin{1};
end
trues=trues(mask>0);


preds=X*Y';
preds=preds(mask>0);

if length(varargin)>3
    % this is a way to evaluate whether the reference data is binary or +/-1.
    shift=varargin{5};
    if shift
        trues=(trues+1)/2.;
    end
end

if length(varargin)>5
    % this is a way to evaluate whether the reference data is binary or +/-1.
    binarise=varargin{7};
    if binarise
        preds=preds>0;
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(error_type,'prscore')
    
end


if strcmp(error_type,'binary')
    err=1-sum(preds == trues)/length(preds);
end

if strcmp(error_type,'l2')
    err=norm(preds-trues,2);
end

if strcmp(error_type,'l2_scaled')
    err=norm(preds-trues,2)/(sqrt(length(preds)));
end

if strcmp(error_type,'cluster')
 
    if length(varargin)<3
        disp('no ground truth cluster provided')
        X=NaN;Y=NaN;
    else
           
        if strcmp(val_idx,'row')
            err=cluster_strength(X,X_true);
        end 
        if strcmp(val_idx,'col')
            err=cluster_strength(Y,Y_true);
        end
    end
end

%For rank_1 recovery:

if strcmp(error_type,'recovered')
    [m,k]=size(X);
    if k>1
        disp('only configured for rank_1')
        err=NaN;
    else
        err=1-(((norm(X-X_true)+norm(Y-Y_true)))>0)*1.;
    end
end
if strcmp(error_type,'recovered_thresh')
    [m,k]=size(X);[n,k]=size(Y);bin_thresh=0.1;
    if k>1
        disp('only configured for rank_1')
        err=NaN;
    else
        err=1-(((norm(X-X_true)^2+norm(Y-Y_true)^2)/(m+n))>bin_thresh)*1.;
    end
end

if strcmp(error_type,'hamming')
    [m,k]=size(X_true);[n,k]=size(Y_true); bin_thresh=0.1;
    if k>1
        disp('only configured for rank_1')
        err=NaN;
    else
    err=(sum(abs(X-X_true))+sum(abs(Y-Y_true)))/(n+m);
    end
end


end
