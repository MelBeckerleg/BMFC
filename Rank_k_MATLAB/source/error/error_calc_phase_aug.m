function [err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming,recovered_rank,num_rows_recovered,num_rows_recovered_95,num_rows_recovered_50,num_rows_recovered_80,num_rows_recovered_75] =error_calc_phase_aug(X,Y,prob,trial_idx,method_idx,k_idx,rho_idx,err_train,err_test,row_cluster_strength,col_cluster_strength,err_hamming,recovered_rank,num_rows_recovered,num_rows_recovered_95,num_rows_recovered_50,num_rows_recovered_80,num_rows_recovered_75)
    [m,k] =size(X);[n,k]= size(Y);num_train = sum(sum(prob.mask));num_test = m*n-num_train;
    preds = X*Y';trues=prob.X_true*prob.Y_true';mask=prob.mask;
    %err_train(trial_idx,method_idx,k_idx,rho_idx) = 1/(2*num_train)*binary_rank_error('l2',[],prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true)^2;    
    %err_test(trial_idx,method_idx,k_idx,rho_idx) = 1/(2*num_test)*binary_rank_error('l2',[],1-prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);%^2;
    err_train(trial_idx,method_idx,k_idx,rho_idx) = nnz((preds(find(mask))>0)*1. - trues(find(mask)))/nnz(mask);%binary_rank_error('binary',[],prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);%^2;    
    err_test(trial_idx,method_idx,k_idx,rho_idx) = nnz((preds(find(1-mask))>0)*1. - trues(find(1-mask)))/nnz(1-mask);%binary_rank_error('binary',[],1-prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);%^2;    
    [pcerr,pcerr_col] = per_cluster_error(X,Y,prob.A_true) ; 
    row_cluster_strength(trial_idx,method_idx,k_idx,rho_idx) = pcerr;
    col_cluster_strength(trial_idx,method_idx,k_idx,rho_idx)= pcerr_col;
    %err_hamming(trial_idx,method_idx,k_idx,rho_idx)=binary_rank_error('hamming','all', prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);
    err_hamming(trial_idx,method_idx,k_idx,rho_idx)=nnz((X*Y'>0)*1. - prob.A_true)/(n*m);%binary_rank_error('binary','all', ones(size(prob.mask)),X,Y,prob.A_true,prob.X_true,prob.Y_true); %I COULD USE THIS TO FIND l2 ERROR
    recovered_rank(trial_idx,method_idx,k_idx,rho_idx) = min(length(find(sum(Y))),length(find(sum(X))));
    
    num_rows_recovered(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0);
    num_rows_recovered_95(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.05);
    num_rows_recovered_50(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.5);
    num_rows_recovered_80(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.2);
    num_rows_recovered_75(trial_idx,method_idx,k_idx,rho_idx) = row_recovered(X*Y',prob.X_true*prob.Y_true',0.25);
end

function [proportion] = row_recovered(A,B,thresh)
    [m,n] = size(A);
    prop_different_by_row = 1/n*sum(abs(A-B),2);
    proportion = 1/m*sum(prop_different_by_row <= thresh);
end