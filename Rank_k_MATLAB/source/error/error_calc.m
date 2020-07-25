function [err_train,err_test,err_row_cluster,err_col_cluster,err_hamming] =error_calc(X,Y,prob,trial_idx,method_idx,err_train,err_test,err_row_cluster,err_col_cluster,err_hamming)
    [m,k]=size(X);[n,k]=size(Y);num_train=sum(sum(prob.mask));num_test=m*n-num_train;
    err_train(trial_idx,method_idx)=1/num_train^2*binary_rank_error('l2','test',prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true)^2;
    err_test(trial_idx,method_idx)=1/num_test^2*binary_rank_error('l2','test',prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true)^2;
    err_row_cluster(trial_idx,method_idx)=binary_rank_error('cluster','row',prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);
    err_col_cluster(trial_idx,method_idx)=binary_rank_error('cluster','col',prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);
    err_hamming(trial_idx,method_idx)=binary_rank_error('hamming','all', prob.mask,X,Y,prob.A_true,prob.X_true,prob.Y_true);
end