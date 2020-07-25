function [Xfinal,Yfinal] = splitting(W_omega,opts)


method=opts.rank_1_method;it=opts.it;k=opts.k;pf=opts.pf;

[Xfinal,Yfinal]=split_op(W_omega,opts);
kval=min(length(find(sum(Xfinal))), length(find(sum(Yfinal))));
kval_final=max(kval,1);
Xfinal=Xfinal(:,1:kval_final);Yfinal=Yfinal(:,1:kval_final);


%allows for columnwise splitting (no improvement observed)
if pf
    kval_final_extra=kval_final;
    for iter=1:kval_final
        [sub_Y,sub_X]=split_op(W_omega(find(Xfinal(:,1)),:)',opts);
        [~,new_k]=size(sub_X);
        
        if kval_final_extra+new_k<k
            extra_X=zeros(m,new_k);
            extra_X(:,Xfinal(:,1))=sub_X;
            extra_Y=sub_Y;
            Xfinal=[Xfinal(:,2:end) extra_X];
            Y_final=[Xfinal(:,2:end) extra_Y];
            kval_final_extra=kval_final_extra+1;
        end
    
    end
    opts.k=kval_final_extra;
end

end
%need to write a rank_k function that can take splitting and then do the
%post processing step on it. 

%function [y,x] TBMC_post(W_omega(find(X(:,split_val),:)')  generates
