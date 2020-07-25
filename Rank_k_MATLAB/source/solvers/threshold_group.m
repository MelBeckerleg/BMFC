function [X] = threshold_group(Mat,k,opts)
%Function for grouping rows within  a tolerance threshold together. 
%Inputs: 
%   Mat (can be either binary or pm1)
%   k - desired rank
%   opts.thresh - Threshold for grouping rows together
%
%Outputs:
%   X - row membership matrix
T=opts.thresh;
[m,n] = size(Mat);

H=squareform(pdist(Mat)); %Not the most memory efficient way to do this

clustered_rows = zeros(m,1);

for k_idx=1:k-1
    unclustered_rows = find(clustered_rows-1);
    if(unclustered_rows)
        selected_row_idx = unclustered_rows(random('Unid',length(unclustered_rows)));
        within_threshold = find(H(selected_row_idx,:)<T);
        X(within_threshold,k_idx) = 1;
        clustered_rows(within_threshold)=1;
    end
end

X(find(clustered_rows-1),k)=1;


end