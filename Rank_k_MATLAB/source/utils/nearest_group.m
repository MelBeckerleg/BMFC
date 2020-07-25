function [X] = nearest_group(Mat,k,opts)
%Function for grouping rows together by iteratively selecting a row and ...
%grouping it with its nearest rows. 
%Inputs: 
%   Mat (can be either binary or pm1)
%   k - desired rank
%
%Outputs:
%   X - row membership matrix

[m,n]=size(Mat);


%allow for the option of observed hamming distance:
observed_hamming=0;
if observed_hamming
    nan_mask = 1-abs(Mat);
    A_copy = Mat;
    A_copy(nan_mask) = Nan;
    H = squareform(pdist(A_copy,@nanhamdist));
else
    H = squareform(pdist(Mat,'hamming'));
end
X = zeros(m,k);


H = squareform(pdist(Mat));
labels = zeros(m,1);
label = 1;

hvals=ones(m,1);
    for j=1:k
        Hsub = H(find(hvals),:);
        [val,idx] = sort(Hsub(j,:));
        idx = idx(1:m/k);
        C = idx;
        for i = 1:length(C)
            iidx = C(i);
            X(iidx,j) = 1;
            hvals(iidx) = 0;
        end
    end
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
