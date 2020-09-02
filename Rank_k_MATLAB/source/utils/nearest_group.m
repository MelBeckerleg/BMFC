function [X_final] = nearest_group(Mat,k,opts)
%Function for grouping rows together by iteratively selecting a row and ...
%grouping it with its nearest rows. 
%Inputs: 
%   Mat (can be either binary or pm1)
%   k - desired rank
%
%Outputs:
%   X - row membership matrix

[m_orig,n_orig]=size(Mat);

% remove almost empty rows
thresh = 0.09; %this should be less than tau but greater than epsilon
empty = find(sum(Mat>0,2) < n_orig*thresh);
not_empty = find(sum(Mat>0,2) >= n_orig*thresh);
Mat = Mat(not_empty,:);
[m,~]=size(Mat);

%allow for the option of observed hamming distance:
observed_hamming=0;
if observed_hamming
    nan_mask = (1-abs(Mat))>0;
    A_copy = Mat;
    A_copy(nan_mask) = nan;
    H = squareform(pdist(A_copy,@nanhamdist));
else
    H = squareform(pdist(Mat,'hamming'));
end
X = zeros(m,k);


%H = squareform(pdist(Mat));
%first =zeros(1);

hvals=ones(m,1);
    for j=1:k
        %first unclassified row:
        unclassified = find(hvals);
        first = unclassified(1);
        distances = H(first,:);
        unclassified_distances = distances(:,[unclassified]);
        %what happens if j>the number of unclassified rows? am I assuming
        %that k<m/k? how do I handle it if not?
        [val,idx] = sort(unclassified_distances);  
        idx = idx(1:floor(m/k)); %do you need to compensate for rounding down here? YES
        C = unclassified(idx);
        for i = 1:length(C)
            iidx = C(i);
            X(iidx,j) = 1;
            hvals(iidx) = 0;
        end
    end
    
X_final = zeros(m_orig,k);
X_final(not_empty,:) = X;
end


% for i=1:m
%     if labels(i)==0
%         [val,idx]=sort(H(i,:));
%         idx=idx(1:ceil(m/k/2));
%         if sum(labels(idx))
%             existing_labels=nonzeros(labels(idx));
%             labels(idx)=existing_labels(1);
%         else
%             labels(idx)=label;
%             label=label+1;
%         end
%     end
% end
% unique_labels=unique(labels);
% recovered_rank=min(length(unique_labels),k);
% for j=1:recovered_rank
%     X(find(labels==unique_labels(j)),j)=1;
% end
