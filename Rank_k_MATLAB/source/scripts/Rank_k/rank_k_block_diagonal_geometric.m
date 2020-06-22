addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/models'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/solvers'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/error'))
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/source/utils'))
clear all; clc;
%opts.epsilon=0.0;opts.k=4;  %opts.rho;opts.a;
opts.epsilon=0.0;
%proportion of allowed incorrect entries - this is different to proportion
%of incorrectly identified rows as I haven't implemented the relative
%cluster score metric ye
thresh=0.03;
%plots for rho/k with tau<1/sqrt(2)

% maximm k occurs when a=1/sqrt(2), 
%you have tau=(1-a)/(1-a^(k));
%and you want the max k such that a^k*tau*n=a^k(1-a)/(1-a^(k))n is greater
%than 1.
% for k_max=1:ceil(sqrt(n))
% val=a^k*(1-a)/(1-a^k)*n;
% if val<1
% k_max=k_max-1;
% break
% end
% end


k_max=4;%12
for m=[50,100,200,300]
    n=m;
opts.m=m;opts.n=m;
a_vals=[0:0.1:0.7 0.71:0.01:0.8 0.9:0.1:1];
num_trials=100;rho_vals=0:0.1:0.9;k_vals=1:k_max;

error=zeros(num_trials,length(a_vals),length(rho_vals),length(k_vals));
error_bin=zeros(num_trials,length(a_vals),length(rho_vals),length(k_vals));
error_bin_thresh=zeros(num_trials,length(a_vals),length(rho_vals),length(k_vals));
error_bin_thresh_rc=zeros(num_trials,length(a_vals),length(rho_vals),length(k_vals));
error_bin_thresh_rc2=zeros(num_trials,length(a_vals),length(rho_vals),length(k_vals));


for trial_idx=1:num_trials
    for rho_idx=1:length(rho_vals)
        opts.rho=rho_vals(rho_idx);
        for a_idx=1:length(a_vals) 
            opts.a=a_vals(a_idx);
            for k_idx=k_max%1:length(k_vals)
            %opts.k=k_vals(k_idx);
            opts.k=k_idx;
            [X_true,Y_true,A_true,W_omega] = generate_block_diagonal_clusters_geometric(opts);
            solve_opts.rank_1_method='lp';solve_opts.it=0;solve_opts.k=opts.k;solve_opts.pf=0;
            [Xfinal,Yfinal] = splitting(W_omega,solve_opts);


            error(trial_idx,a_idx,rho_idx,k_idx)=norm(A_true-Xfinal*Yfinal');
            error_bin(trial_idx,a_idx,rho_idx,k_idx)=1*(norm(A_true-Xfinal*Yfinal')>0);
            error_bin_thresh(trial_idx,a_idx,rho_idx,k_idx)=1*(norm(A_true-Xfinal*Yfinal')>(thresh*sqrt(m)*sqrt(n)));
            %error_bin_thresh_rc(trial_idx,a_idx,rho_idx)=1*((norm(X_true-Xfinal)+norm(Y_true-Yfinal))>(thresh*sqrt(m)*sqrt(n)*opts.k));
            %if opts.k<5
             %   cm=cluster_metric();
              %  error_bin_thresh_rc2(trial_idx,a_idx,rho_idx)=1*((cm.cluster_error_small_r(X_true,Xfinal,opts.k)+cm.cluster_error_small_r(Y_true,Yfinal,opts.k))>(thresh*sqrt(m)*sqrt(n)*opts.k));
            %end
            %save(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotau.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh')
            save(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotau_highres_100trials.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh')
            %save(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotauk.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh')
            clear W_omega W_omega_new A_true
            end
        end
    end
end
end

% @(k)(1-2*0.7.^2+0.7.^(2*k))./(3+2*0.7^2-0.7.^(2*k))
 
 %t1= @(a,k) (1-a)./(1-a^k);
 %mu1=@(n,rho,a) t1(a,k).^2.*n.^2.*rho./2;
 %mu2=@(n,rho,a) t1(a,k).^2.*n.^2.*rho./2.*((1-a.^k)./(1-a)-1);
 %probA=@(delta,n,a,rho,k) 1-exp(-delta.^2.*mu1(n,rho,a)./(2+delta));
 %probB=@(delta,n,a,rho,k) 1-exp(-delta.^2.*mu2(n,rho,a)./(2));
 %prob=@(delta,n,a,rho,k)  probA(delta,n,a,rho,k).*probB(delta,n,a,rho,k);

 %delta_fun= @(a,k)  1-2.*(1-a)./(1-a.^(2.*k));
 %k_d=delta_fun(a,k_vals);
if 0
    load(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotauk.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh') 
    figure();[t,a,r,c]=size(error_bin);surf(reshape(mean(error_bin),r,c)');view([0,90])
     my_mat_min=zeros(length(rho_vals),length(k_vals));my_mat_max=my_mat_min;
     for rr=1:length(rho_vals)
         for kk=1:length(k_vals)
    my_mat_min(rr,kk)=min_success_prob(300,0.7,rho_vals(rr),k_vals(kk));
    my_mat_max(rr,kk)=max_success_prob(300,0.7,rho_vals(rr),k_vals(kk));

    %my_mat_min(rr,kk)=min_success_prob_independent(300,0.7,rho_vals(rr),k_vals(kk));
    %my_mat_max(rr,kk)=max_success_prob_independent(300,0.7,rho_vals(rr),k_vals(kk));
         end
     end
     figure()
     surf(my_mat_max);view([0,90]);title('max');shg
     colorbar
     figure()
     surf(my_mat_min);view([0,90]);title('min');shg
     colorbar
 
 end
if 0
    plot_idx=1;
    s1=figure();s2=figure();s3=figure();
    figure_labels={'s1','s2','s3'};
    for m=[50,100,250]
    %load(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotau.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh')
    load(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotau_highres.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh');error=error(:,:,:,end);error_bin=error_bin(:,:,:,end);error_bin_thresh=error_bin_thresh(:,:,:,end);
    [t,r,c]=size(error);
    set(0,'defaulttextinterpreter','latex')%,'defaulttextfontsize',10,'DefaultAxesFontSize',10)
    
    for error_type={'error','error_bin','error_bin_thresh'}
    h=surf(reshape(mean(eval(error_type{1})),r,c)'); view([0,90]);set(h,'EdgeColor','None')
    
    xticklabels(a_vals);xlabel('a');yticklabels(rho_vals);ylabel('$\rho$')
    set(gca,'FontSize',20);xlim([1,length(a_vals)]);ylim([1,length(rho_vals)]);
    %saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/geometric_m=%s_%s',string(m),string(error_type)),'epsc')
    end
    end
    
    
    s1=figure();s2=figure();s3=figure();figure_labels={'s1','s2','s3'};
    plot_idx=1;
    for m=[50,100,200,300]%[50,100,125,150,200,250,300]   
    load(sprintf('/home/user/Documents/Mel/Ethera/Results/geometric_clusters_m%s_differentarhotau_highres_100trials.mat',string(m)),'error','error_bin','error_bin_thresh','error_bin_thresh');
    error=error(:,:,:,end);error_bin=error_bin(:,:,:,end);error_bin_thresh=error_bin_thresh(:,:,:,end);
    [t,r,c]=size(error);
    error_val=1;
    for error_type={'error_bin_thresh','error_bin'}%{'error','error_bin','error_bin_thresh'}
    figure(eval(figure_labels{error_val}))
    
    subplot(1,4,plot_idx)
    val=reshape(mean(eval(error_type{1})),r,c)'
    idx=[vec(repmat(1:7,10,1))' 8:18 vec(repmat(19:20,10,1))'];
    %idx=1:length(a_vals)
    h=surf(val(:,idx)); view([0,90]);set(h,'EdgeColor','None')
    title(sprintf('m=%s',string(m)),'Interpreter','Latex');
    x_vals=a_vals(idx);
    xticks(1:20:length(x_vals));yticks(1:2:length(rho_vals));
    xticklabels(x_vals(1:20:length(x_vals)));xlabel('a','Interpreter','Latex');yticklabels(rho_vals(1:2:length(rho_vals)));ylabel('$\rho$','Interpreter','Latex')
    line_point=find(a_vals(idx)>=0.72);line_point=line_point(1);
    hold on; line([line_point+1;line_point+1],[0; 10],[1,1],'LineWidth',3,'LineStyle','--','Color','k')
    set(gca,'FontSize',10);xlim([1,length(x_vals)]);ylim([1,length(rho_vals)]);
    saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/geometric_%s_highres_100trials',string(error_type)),'epsc')
    
    
    
    
    error_val=error_val+1;
    end 
    
    plot_idx=plot_idx+1;
    
    end
end
%% Investigate interactive plotting for different values of m

%     f=figure()
%     ax=axes('Parent',f,'position',[0.13 0.39  0.77 0.54])
%     h=surf(zeros(11,10));
%     xticklabels(a_vals);xlabel('a');yticklabels(rho_vals);ylabel('$\rho$')
%     set(gca,'FontSize',20);xlim([1,length(a_vals)]);ylim([1,length(rho_vals)]);
% 
%     b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
%               'value',M, 'min',50, 'max',150);
%     
% 
%     bgcolor = f.Color;
%     bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
%                 'String','1','BackgroundColor',bgcolor);
%     bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
%                 'String','150','BackgroundColor',bgcolor);
%     bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
%                 'String','Number of rows','BackgroundColor',bgcolor);
% 
%     
%     b.Callback = @(m) updateSystem(h,load_m_error(m)); 
%     
%     
%     
% end
% 
% function load_m_error(m)
% 
% end
