function [Xfinal,Yfinal]=split_op(W_omega,opts)

method=opts.rank_1_method;it=opts.it;k=opts.k;
[m,n]=size(W_omega);

max_store_val=3*k;%This stops too many partitions being stored. It would be interesting to
%adapt the tolerance if too many partitions were being generated....



Xstore=zeros(m,2*k); %store the partitions to try
Xfinal=zeros(m,2*k); %when a partition block cannot be subdivided, it is added to Xfinal, Yfinal
Yfinal=zeros(n,k);

splits=0;
Xstore(:,1)=ones(m,1);
store_val=1;kval=1;
while kval<=k+1 && store_val>=1 && sum(sum(Xstore)) ;
    %disp('loop')
    [x,y]=rank_1_solve(W_omega(Xstore(:,1)>0,:),method);
    %AVG DOESN'T WORK IF x has dimension 1! possibly the same for partition
    if it
        [x,y]=iterative_update(x,y,W_omega(Xstore(:,1)>0,:));
    end
    %ITERATIVE UPDATING WITH the wrong dimensions doesn't work!
    %check binary
    %if this has changed 
    %% Section 
    if sum(1-x) %if there are 
        if sum(x) && sum(y)
            if store_val+1>2*k
                disp('more partitions to be check than allocated')
            end
            
            if store_val+2<max_store_val
            %store both sides of the partition to be considered further
            positive_side=Xstore(:,1);
            positive_side(positive_side>0)=x;
            negative_side=Xstore(:,1);
            negative_side(negative_side>0)=1-x;
            Xstore(:,store_val+1)=positive_side;
            Xstore(:,store_val+2)=negative_side;   
            %future patterns need to not overlap with existing stored
            %patterns (note that this will be dropped by 1 later)
            store_val=store_val+2;
            else
                disp('too many partitions to be checked')
            end
        end

        
    else
        if sum(y)%if no further partitioning suggested, and a 
            %nonzero tile is recovered, save the tile
            %disp('storing')
            Xfinal(:,kval)=Xstore(:,1);
            Yfinal(:,kval)=y;       
            kval=kval+1;              
        end
    end
    
    
    %drop the pattern that you have just explored, future patterns need to
    %be stored one column earlier
    
    Xstore=Xstore(:,2:end);
    store_val=store_val-1;
    
        
end