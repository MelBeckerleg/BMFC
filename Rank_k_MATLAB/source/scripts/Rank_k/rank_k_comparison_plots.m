addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/'))

%Combine test for different m
test_params='num_trials=30;epsilon=0.03;tau=0.3;rho=0.7;k=@(x) log(x);...'
                '...num_methods=14;m_vals = [50,100,150,200];';

                
err_test_results = zeros(num_trials,num_methods,length(m_vals));
err_train_results = zeros(num_trials,num_methods,length(m_vals));
comp_time_results = zeros(num_trials,num_methods,length(m_vals));
col_cluster_strength_results = zeros(num_trials,num_methods,length(m_vals));
row_cluster_strength_results = zeros(num_trials,num_methods,length(m_vals));



%load for different mvalues

iter=1;

for m = m_vals
    k=floor(log(m));
    fn=sprintf('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_%s_k_%s_rho_%s_tau_%s_eps_%s.mat',string(m),string(k),string(100*rho),string(100*tau),string(100*epsilon));
    load(fn);
    
%     
%     err_test_results(iter,:) = (err_test);
%     err_train_results(iter,:) = mean(err_train);
%     comp_time_results(iter,:) = mean(comp_time);
%     col_cluster_strength_results(iter,:) = mean(col_cluster_strength);
%     row_cluster_strength_results(iter,:) = mean(row_cluster_strength);
%     err_test_std(iter,:) = std(err_test);
%     err_train_std(iter,:) = std(err_train);
%     comp_time_std(iter,:) = std(comp_time);
%     col_cluster_strength_std(iter,:) = std(col_cluster_strength);
%     row_cluster_strength_std(iter,:) = std(row_cluster_strength);

        err_test_results(:,:,iter) = err_test;
        err_train_results(:,:,iter) = (err_train);
        comp_time_results(:,:,iter) = (comp_time);
        col_cluster_strength_results(:,:,iter) = (col_cluster_strength);
        row_cluster_strength_results(:,:,iter) = (row_cluster_strength);


    iter=iter+1;
end



method_labels={};iter=1;
for new_name = names
    
    method_labels{iter}=char(join(new_name{1}));
    iter = iter+1;
end


fn=sprintf('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_all_m_.mat')
save(fn,'m_vals','names','err_test_results', 'err_train_results','comp_time_results','col_cluster_strength_results','row_cluster_strength_results','method_labels','params');

if 0:
    %This is a very fixed way of setting the formatting, be careful if you change any of the methods.
    markers = {'*','<','+','x','o','p','d','h','v','^','>','.','s','none'};
    linestyles={':','-',':','-',':','-',':','-','--','-.','-.','-.','-.','-'};
    colors={'b','b','g','g','g','g','c','c','m','k','k','k','k','r'};
    method_idx=[1:num_methods];
    mmarkers=containers.Map(method_labels,markers);
    mlinestyles=containers.Map(method_labels,linestyles);
    mcolors=containers.Map(method_labels,colors);
    midx=containers.Map(method_labels,method_idx);
    
    %CHECK THIS STEP FOR SURE!
    names = {'TBMC','TBMC-AM','Split-avg','Split-avg-AM','Split-pivot','Split-pivot-AM','Split-bnmf','Split-bnmf-AM','Neighbour','NuclearNorm','Convex','1bit','!LRS','BNMF'};
    mnames=containers.Map(method_labels,names);
    
    
    
    %Now let's set the maps for the different error metrics
    results = {'err_test_results', 'err_train_results','comp_time_results','col_cluster_strength_results','row_cluster_strength_results'};
    result_titles = {'Test Error', 'Training Error', 'Computational Time', 'Column Recovery','Row Recovery'};
    result_names_y = {'Proportion of entries missclassified','Proportion of entries missclassified','Time (s)', 'Strength','Strenth'};
    mresult_names_y = containers.Map(results,result_names_y);
    mresult_titles = containers.Map(results,result_titles);
end

if 0:
    tag='mid_july';
    save(sprintf('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_all_m_method_configs_%s_%s.mat',tag),'params');
    methods=method_labels[1:4,6:end];%choose a subset if you like.
    for result=results
        close('all');figure('units','normalized','outerposition',[0 0 1 1]);
        plot_val=eval(result{1});
        plot_val_mean=squeeze(mean(plot_val));plot_val_std=squeeze(std(plot_val));
        close('all')
        for method = methods%alternatively you can select the methods that you want to plot
        hold on
        errorbar(m_vals,plot_val_mean(midx(method{1}),:),plot_val_std(midx(method{1}),:),'LineWidth',2,'Marker',mmarkers(method{1}),'LineStyle',mlinestyles(method{1}),'MarkerSize',10,'Color',mcolors(method{1}))
        xlabel('m');title(mresult_titles(result{1}));ylabel(mresult_names_y(result{1}));
        shg
        end
        %check this legend, you're using the mapped labels
        legend(values(mnames,methods),'Location','bestoutside')
        saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_all_m_method_%s_%s.epsc',result{1},tag))
    end
    
end