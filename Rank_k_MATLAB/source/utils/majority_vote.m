function [Y] = majority_vote(Mat,X)
% takes in a clustering and matrix in {1,-1,0}, where 0 denotes a missing 
% entry, and outputs Y, the row footprints for each cluster
[m,n]=size(Mat);
[~,k]=size(X);
Y=zeros(n,k);

for j=1:k
    idx=find(X(:,j));
    mu=(sum(Mat(idx,:))>0)*1;
    Y(:,j)=mu;
end


end