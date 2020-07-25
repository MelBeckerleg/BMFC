function [X,Y] = spectral(M,opts)
    [m,n]=size(M);k=opts.k;
    %erasure probability
    epsilon=1.-nnz(M)/(1.*m*n);
    %threshold for distance from user.  /!\
    thresh=sqrt(1-epsilon)*k*log(n)*12;
    %assign subsets
    %probability of assignment
    delta=(1-epsilon)/4.;
    omega=find(M);
    rho_approx=1.*length(omega)/m/n;
    assign=rand(length(omega),1);
    omega_one=find(assign<(1/2.));
    omega_two=find((assign>(1/2.-delta)).*(assign<(1.-delta)));

    R_one=zeros(m,n);R_two=zeros(m,n);
    R_one(omega(omega_one))=M(omega(omega_one));
    R_two(omega(omega_two))=M(omega(omega_two));

    [U,S,V]=svd(R_one);
    [mr,nr]=size(S);
    if min(mr,nr)>k
        PR_one=1.*(U(:,1:k)*S(1:k,1:k)*V(:,1:k)');
    else
        PR_one=U*S*V';
    end

    [U,S,V]=svd(R_two);
    [mr,nr]=size(S);
    if min(mr,nr)>k
        PR_two=1.*(U(:,1:k)*S(1:k,1:k)*V(:,1:k)');
    else
        PR_two=U*S*V';
    end


    to_cluster=ones(m,1);
    row_clusters=zeros(m,1);
    
    iter=0;
    while iter<k && sum(to_cluster)
        options=find(to_cluster);
        row=options(randi(length(find(to_cluster))));
        chosen=PR_one(row,:);
        indices=zeros(m,1);
        for ii=1:m
            indices(ii)=1.*(norm(PR_one(ii,:)-chosen,2)^2<thresh);
        end
	indices=find(indices);
    to_cluster(indices)=0;
	row_clusters(indices)=iter+1;
    iter=iter+1;
    end
    k_row=iter;
    %arbitrary assignment of the remaining rows
    row_clusters(find(to_cluster))=randi(k_row,length(find(to_cluster)),1);

    to_cluster=ones(n,1);
    col_clusters=zeros(n,1);

    iter=0;
    while iter<k && sum(to_cluster)
        options=find(to_cluster);
        col=options(randi(length(find(to_cluster))));
        chosen=PR_one(:,col);
        indices=zeros(n,1);
        for ii=1:n
            indices(ii)=1.*(norm(PR_one(:,ii)-chosen,2)^2<thresh);
        end
	indices=find(indices);
    to_cluster(indices)=0;
	col_clusters(indices)=iter+1;
    iter=iter+1;
    end
    
    k_col=iter;


    %arbitrary assignment of the remaining cols
    col_clusters(find(to_cluster))=randi(k_col,length(find(to_cluster)),1);

    %block assign values

    
    %recluster
    %find the mean of the rows 
    row_centres=zeros(k_row,n);
    for ii=1:k
        row_centres(ii,:)=mean(PR_two(find(row_clusters==ii),:));
    end

    %create distance
    distance=zeros(m,k_row);
    for ii=1:m
    for jj=1:k_row
        distance(ii,jj)=sum(abs(row_centres(jj,:)-PR_two(ii,:)));
    end
    end

    min_vals = find((distance-min(distance,[],2))==0 );
    X=zeros(m,k_row);
    for ii=1:k_row
        X(:,ii)=1.*(min_vals==ii);
    end

    % find the mean of the columns
%     col_centres=zeros(k,m);
%     for ii=1:k
%         col_centres(ii,:)=mean(PR_two(:,find(col_clusters==ii))');
%     end
% 
% %     %create distance
% %     distance=zeros(n,k);
%     for ii=1:n
%     for jj=1:k
%         distance(ii,jj)=norm(row_centres(jj,:)-PR_two(:,ii)');
%     end
%     end
% 
%     min_vals = find(distance-min(distance,[],2)==0 );
%     Y=zeros(n,k);
%     for ii=1:k
%         Y(:,ii)=1.*(min_vals==ii);
%     end
    
    %assign block values
        block=zeros(k_row,k_col);
    for ii=1:k_row
    for jj=1:k_col
        row_idx=find(row_clusters==ii);
        col_idx=find(col_clusters==jj);
        block(ii,jj)=(sum(sum(PR_two(row_idx,col_idx)))>0)*1.;
    end
    end
    
  
    Y=zeros(n,k_row);
    for ii=1:k_row
        for pos_col_cluster=find(block(ii,:))
            Y(:,ii)=Y(:,ii)+1.*(col_clusters==pos_col_cluster);
        end
    end
        
    
    
    
end
    







