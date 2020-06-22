close all;clear;
fn_home='/home/user/Documents/Mel/Ethera/Results/'
save_label='two_block_approx';
%fn_home_im='/home/user/Documents/Mel/Ethera/FirstYear/Images/';fn_home='/home/user/Documents/Mel/Ethera/SecondYear/EvilPlans/';
%fn_home_im='/home/beckerleg/Ethera/FirstYear/Images/';fn_home='/home/beckerleg/Ethera/FirstYear/Clusters/SecondYear/PhaseWork/';
bin_thresh=0.05; %Ideally this would be 0
m = 10;
n = 10;
num_trials=3;
step_val=0.1;tau_vals=0.1:step_val:1;tv=length(tau_vals);
rho_vals=0:0.1:1;rv=length(rho_vals);eps=0;epsilon=eps;
objective=zeros(tv,tv,rv);error=zeros(tv,tv,rv);error_bin=zeros(tv,tv,rv);error_bin_thresh=zeros(tv,tv,rv);

%for tau1=tau_vals
l=1;
for rho=rho_vals
i=1;j=1;
    for tau2=tau_vals
        tau1=tau2;
        for rand_trials=1:num_trials
            scaled_n=n/(1+tau1);
            scaled_m=m/(1+tau1);
            x_true=zeros(m,1);y_true=zeros(n,1);
            x_true(1:max(floor(scaled_m),1))=1;y_true(1:max(floor(scaled_n),1))=1; 
            x_true2=zeros(m,1);y_true2=zeros(n,1);
            x_true2(min(floor(scaled_m)+1,m):end)=1;y_true2(min(floor(scaled_n)+1,n):end)=1; 
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
            fn=sprintf('%stwo_block_error_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'error')
            fn=sprintf('%stwo_block_errorbin_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'error_bin')
            fn=sprintf('%stwo_block_errorbinthresh_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'error_bin_thresh')
            fn=sprintf('%stwo_block_errorobjective_epsilon%s',fn_home,string(100*epsilon));
            save(fn,'objective')
        end
        error(i,j,l)=error(i,j,l)/num_trials;
        error_bin(i,j,l)=error_bin(i,j,l)/num_trials;
        error_bin_thresh(i,j,l)=error_bin_thresh(i,j,l)/num_trials;
        j=j+1;
    end
    %i=i+1; j=1;
%end
i=1;j=1;l=l+1;
	end



%%%%%%%%%%%%Plots




 %% Plots on my computer
tau_vals=0.1:step_val:1;
rho_vals=0:0.1:1;
fn_home_im='/home/user/Documents/Mel/Ethera/FirstYear/Images/'
fn_home='/home/user/Documents/Mel/Ethera/Results/'%folder temp for ntrials=3
set(0,'defaulttextinterpreter','latex')
for epsilon=eps
    
    fn=sprintf('%stwo_block_error_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'error')
    fn=sprintf('%stwo_block_errorbin_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'error_bin')
    fn=sprintf('%stwo_block_errorbinthresh_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'error_bin_thresh')
    fn=sprintf('%stwo_block_errorobjective_epsilon%s',fn_home,string(100*epsilon));
    load(fn,'objective')
  
    
    for error_type={'hamming','recovered','objective','threshrecovered'}
        if strcmp(error_type,'recovered')
            plot_var=error_bin;
            label_for_y_axis='Proportion recovered';

        elseif strcmp(error_type,'hamming')
            plot_var=error;
            label_for_y_axis='Distance from optimal';

        elseif strcmp(error_type,'objective')
            plot_var=objective;
            label_for_y_axis='Objective';
        elseif strcmp(error_type,'threshrecovered')
            plot_var=error_bin_thresh;
            label_for_y_axis='Proportion recovered';
        end

        %for fixed rho=0 (no missing data), plot the recovery space
       
        
        if 1
            figure();hold on;
            for r_idx=1:rv
                markers={'^','o','+','>','<','*','s','d','p','','.'};
                vals=plot_var(1,:,r_idx);
                plot(tau_vals(2:end),vals(2:end),sprintf('-%s',markers{r_idx}),'DisplayName',sprintf('%s=%0g\n',string('\rho'),round(rho_vals(r_idx),2,'significant')),'LineWidth',2,'MarkerSize',10)
                %names{r_idx}=sprintf('%s=%0g\n',string('\rho'),round(rho_vals(r_idx),2,'significant'));
                hold on;
                xlabel('$\tau$')
                ylabel(label_for_y_axis)
                fn=sprintf('%stwo_block_tau1tau2recoveryproportionepsilon%s%soct',fn_home_im,string(100*epsilon),string(error_type));
                legend('AutoUpdate','on')
                
                saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');
            end
            hold off; 
            %legend(gca,names)
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
