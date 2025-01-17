close all;clear;
fn_home='/home/user/Documents/Mel/Ethera/Results/'
save_label='row_overlap_tmp';
%fn_home_im='/home/user/Documents/Mel/Ethera/FirstYear/Images/';fn_home='/home/user/Documents/Mel/Ethera/SecondYear/EvilPlans/';
%fn_home_im='/home/beckerleg/Ethera/FirstYear/Images/';fn_home='/home/beckerleg/Ethera/FirstYear/Clusters/SecondYear/PhaseWork/';
bin_thresh=0.05; %Ideally this would be 0
m = 30;
n = 30;
num_trials=5;
step_val=0.05;tau_vals=0.1:step_val:1;tv=length(tau_vals);
rho_vals=0:0.1:1;rv=length(rho_vals);rho=0;k=1;eps=0;epsilon=eps;
objective=zeros(tv,tv,rv);error=zeros(tv,tv,rv);error_bin=zeros(tv,tv,rv);error_bin_thresh=zeros(tv,tv,rv);
l=1;
for rho=rho_vals  
i=1;j=1;
for tau1=tau_vals
    for tau2=tau_vals
        for rand_trials=1:num_trials
            %we need to scale so that the total sum is n. this means  
            %find N such that N(1+tau-eta*tau)=n
            %so N=n/(1+tau-eta*tau);
            scaled_n=n/(1+tau1-tau1*tau2);
            scaled_m=m/(1+tau1);
            x_true=zeros(m,1);y_true=zeros(n,1);
            x_true(1:max(floor(scaled_m),1))=1;y_true(1:max(floor(scaled_n),1))=1; 
            x_true2=zeros(m,1);y_true2=zeros(n,1);
            x_true2(min(floor(scaled_m)+1,m):end)=1;y_true2(min(max(ceil((1-tau1*tau2)*scaled_n),1),n):end)=1; 
            X=x_true*y_true'+x_true2*y_true2';
            X = X + (1-2*X).*(rand(m,n)<epsilon); 
            W = 2*X-1;
            W_omega=( W+(1-2*X).*(rand(m,n)<rho));


            numzs=length(find(W_omega(:)<0));
            opt=[zeros(numzs,1);x_true;y_true];
            %need to set z=1 wherever (5*x_true*y_true'-W_omega)==6; 
            %and identify the index of z this corresponds to 
            fi=[find((5*x_true*y_true'-W_omega)==6)];
            for flipped_idx=1:length(fi)
                flipped=fi(flipped_idx);
                opt(find(find(W_omega<0)==flipped))=1;
            end

            c = -[-ones(numzs,1); 1/2*sum((W_omega')>0)'; 1/2*sum((W_omega)>0)'];
            %%%this is not the most memory efficient thing to do, but since we need to go up to rho=1 it'll do...
            %%%%we index z such that zk=zij -> zk+1=zi+1j unless i=m.
            %%%%then A is  stacking n copied of eye(m) followed by m copies of n ones. 
            %%%%%1     1  1
            %     1     1 1
            %      1   1   1
            %       1   1  1  
            A = [-eye(n*m) kron(ones(m,1),eye(n)) kron(eye(m),ones(n,1))];
            %check indexing!
            A=A([find(W_omega(:)<0)],:);
            A=A(:,[find(W_omega(:)<0)' (m*n+1:1:n*m+m+n)]);
            b = [ones(numzs,1)];

            options = optimoptions('linprog','Algorithm','dual-simplex');
            out = linprog(c,A,b,[],[],zeros(numzs+m+n,1),ones(numzs+m+n,1),options);
            %Z = reshape(out(1:numzs,1),m,n);
            
            u = out(numzs+1:numzs+m,1);
            v = out(numzs+m+1:end,1);

            %success?
            error_bin(i,j,l)=error_bin(i,j,l)+1-(((norm(u-x_true)+norm(v-y_true)))>0)*1.;
            error_bin_thresh(i,j,l)=error_bin_thresh(i,j,l)+1-(((norm(u-x_true)^2+norm(v-y_true)^2)/(m+n))>bin_thresh)*1.;
            error(i,j,l)=error(i,j,l)+(sum(abs(u-x_true))+sum(abs(v-y_true)))/(n+m);
            c'*(out-opt);
            objective(i,j,l)=objective(i,j,l)+c'*(out-opt);
            fn=sprintf('%s_two_block_error_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'error')
            fn=sprintf('%s_two_block_errorbin_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'error_bin')
            fn=sprintf('%s_two_block_errorbinthresh_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'error_bin_thresh')
            fn=sprintf('%s_two_block_errorobjective_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'objective')
        end
        error(i,j,l)=error(i,j,l)/num_trials;
        error_bin(i,j,l)=error_bin(i,j,l)/num_trials;
        error_bin_thresh(i,j,l)=error_bin_thresh(i,j,l)/num_trials;
        i=i+1;
    end
    i=1; j=j+1;
end
i=1;j=1;l=l+1;
	end



%%%%%%%%%%%%Plots




 %% Plots on my computer
tau_vals=0.1:step_val:1;tau_vals=tau_vals(1:end-1);
rho_vals=0:0.1:1;
fn_home_im='/home/user/Documents/Mel/Ethera/FirstYear/Images/'
fn_home='/home/user/Documents/Mel/Ethera/Results/'%folder temp for ntrials=3
set(0,'defaulttextinterpreter','latex')
for epsilon=eps
    
    fn=sprintf('%s_two_block_error_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'error')
    fn=sprintf('%s_two_block_errorbin_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'error_bin')
    fn=sprintf('%s_two_block_errorbinthresh_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'error_bin_thresh')
    fn=sprintf('%s_two_block_errorobjective_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'objective')
  

    for error_type={'hamming','recovered','objective','threshrecovered'}
        if strcmp(error_type,'recovered')
            plot_var=error_bin;

        elseif strcmp(error_type,'hamming')
            plot_var=error;

        elseif strcmp(error_type,'objective')
            plot_var=objective;
        elseif strcmp(error_type,'threshrecovered')
            plot_var=error_bin_thresh;
        end

        %for fixed rho=0 (no missing data), plot the recovery space
        if 0
            figure()

            imagesc(plot_var(1:end-1,1:end-1,1))
            set(gca,'YDir','normal','Fontsize',20);xlabel('$\tau_1$','Interpreter','Latex');ylabel('$\tau_2$','Interpreter','Latex')
            xticklabels=tau_vals([1,ceil(length(tau_vals)/4),ceil(length(tau_vals)/2),ceil(3*length(tau_vals)/4),length(tau_vals)]) ;
            xticks=linspace(1,size(error,2),numel(xticklabels));
            set(gca,'Xtick', xticks,'XTickLabel',xticklabels,'Ytick',xticks,'YTickLabel',xticklabels,'Fontname','Latin Modern Roman')
            fun1=@(a,e) e*(a-1)./((-1+4*a)*e-a);plotscale=@(truth) (truth-tau_vals(1))*1/step_val+1;
            y=zeros(length(tau_vals),1);i=0;y1=y;y2=y;y3=y;
            for a=tau_vals
                i=i+1;
                y1(i)=plotscale(fun1(a,epsilon));
            end
            hold on; plot(y1,'color','k','linewidth',2);
            val=plotscale(epsilon/(2*(1-epsilon)));
            line([0,numel(tau_vals)],[val,val],'color','k','linewidth',2);line([val,val],[0,numel(tau_vals)],'color','k','linewidth',2);
            xlim([0,numel(tau_vals)]);ylim([0,numel(tau_vals)]);colorbar
            fn=sprintf('%s_two_block_tau1tau2recoveryphaseepsilon%s%soct',fn_home_im,string(100*epsilon),string(error_type));
            saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close


            %Plot how recovery changes with rho for different Tau. 
            markers={'^','o','+','>','<','*','s','d','p','.',''};
            set(0,'defaulttextInterpreter','latex');set(legend,'Interpreter','latex')
            figure()
            for i=1:length(tau_vals)
                hold on 
            legend('AutoUpdate','on')
            plot(rho_vals,reshape(plot_var(i,i,:),1,[]),sprintf('-%s',markers{i}),'DisplayName',sprintf('Tau=%0g\n',round(tau_vals(i),2,'significant')),'LineWidth',2,'MarkerSize',10)
            legend('AutoUpdate','on')
            end
            hold off

            xlabel('$\rho$')
            ylabel(sprintf('Error:%s',string(error_type)))
            fn=sprintf('%s_two_block_tau1tau2recoveryphasemissingepsilon%s%soct',fn_home_im,string(100*epsilon),string(error_type));
            saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close
        end
        
        if 0
            figure();hold on;
            for r_idx=1:rv
                markers={'^','o','+','>','<','*','s','d','p','.',''};
                vals=plot_var(1,:,r_idx);
                plot(tau_vals,vals(2:end),sprintf('-%s',markers{r_idx}),'DisplayName',sprintf('%s=%0g\n',string('\rho'),round(rho_vals(r_idx),2,'significant')),'LineWidth',2,'MarkerSize',10)
                %names{r_idx}=sprintf('%s=%0g\n',string('\rho'),round(rho_vals(r_idx),2,'significant'));
                hold on;
                xlabel('$\tau$')
                ylabel('Proportion solved')
                fn=sprintf('%s_two_block_tau1tau2recoveryproportionepsilon%s%soct',fn_home_im,string(100*epsilon),string(error_type));
                legend('AutoUpdate','on')
                
                saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');
            end
            hold off; 
            %legend(gca,names)
        end
            
        if 1
            
            figure()

            imagesc(plot_var(1:end-1,1:end-1,1))
            set(gca,'YDir','normal','Fontsize',20);xlabel('$\tau$','Interpreter','Latex');ylabel('$\eta$','Interpreter','Latex')
            xticklabels=tau_vals([1,ceil(length(tau_vals)/4),ceil(length(tau_vals)/2),ceil(3*length(tau_vals)/4),length(tau_vals)]) ;
            xticks=linspace(1,size(error,2),numel(xticklabels));
            set(gca,'Xtick', xticks,'XTickLabel',xticklabels,'Ytick',xticks,'YTickLabel',xticklabels,'Fontname','Latin Modern Roman')
            fun2=@(a) 1/a*(1-a^2);
            plotscale=@(truth) (truth-tau_vals(1))*1/step_val-1;
            fun3=@(a) (1-a)/(2*a);
            fun4 = @(a) 2*(1-a)/(1+a);
            %fun3=@(a) a/(a+1);
            
            y=zeros(length(tau_vals),1);i=0;y1=y;y2=y;y3=y;
            for a=tau_vals
                i=i+1;
                y1(i)=plotscale(fun2(a));
                y2(i)=plotscale(fun3(a));
                y3(i)=fun3(a);
                y4(i)=plotscale(fun4(a));
            end
            hold on; plot(y1,'color','k','linewidth',2);
            val=plotscale(epsilon/(2*(1-epsilon)));
            %line([0,numel(tau_vals)],[val,val],'color','k','linewidth',2);line([val,val],[0,numel(tau_vals)],'color','k','linewidth',2);
            ymin=1/sqrt(2)/(1/sqrt(2)+1);
            line([plotscale(1/sqrt(2)),plotscale(1/sqrt(2))],[plotscale(0),numel(tau_vals)],'color','k','linewidth',2,'LineStyle','--')
            coods=find(tau_vals>=1/sqrt(2));cood=coods(1);
            plot(y2(1:cood-3),'color','k','linewidth',2,'LineStyle','--');
            plot(y4(1:cood-3),'color','k','linewidth',2,'LineStyle','--');
            set(groot,'defaultAxesTickLabelInterpreter','latex'); 
            xlim([plotscale(tau_vals(3)),numel(tau_vals)]);ylim([plotscale(tau_vals(3)),numel(tau_vals)]);%colorbar
            fn=sprintf('%s_two_block_tau1tau2recoveryphaseepsilon%s%s%s',fn_home_im,string(100*epsilon),string(error_type),string(save_label));
            saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close

        end 
    end
    
end
% c1=@(tv) 2*tv/(1+2*tv);
% c2=@(tv1,tv2) tv1*tv2/(1+4*tv1*tv2-tv1-tv2);
% min_epsilon=zeros(length(tau_vals),length(tau_vals))
% i=1;j=1;
% for tv1=tau_vals
%     for tv2 = tau_vals
%         min_epsilon(i,j)=min(min(c1(tv1),c1(tv2)),c2(tv1,tv2));
%         i=i+1;end
%     i=1;j=j+1;
% end
% figure()
% imagesc(min_epsilon);fn='/home/user/Documents/Mel/Ethera/FirstYear/Images/minvalepsilonrank1convex'
% ;%surf(min_epsilon)zlabel('$\epsilon$','Interpreter','Latex');fn='/home/user/Documents/Mel/Ethera/FirstYear/Images/minvalepsilonrank1convex_surf'
% set(gca,'YDir','normal');xlabel('$\tau_1$','Interpreter','Latex');ylabel('$\tau_2$','Interpreter','Latex');
% xticklabels=tau_vals([1,ceil(length(tau_vals)/4),ceil(length(tau_vals)/2),ceil(3*length(tau_vals)/4),length(tau_vals)]) ;
% xticks=linspace(1,size(error,2),numel(xticklabels));
% set(gca,'Xtick', xticks,'XTickLabel',xticklabels,'Ytick',xticks,'YTickLabel',xticklabels,'FontName','Times','FontSize',15)
% colorbar()
% 
% saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close
