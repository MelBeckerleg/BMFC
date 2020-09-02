
%% Set-up
addpath(genpath('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/'))
fn='/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_500_tau_30_eps_5_aug_1bit_long.mat';
%fn = '/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_500_tau_30_eps_5_all_aug22_m500_NeighbourNuclear_long.mat'
fn = '/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_500_tau_30_eps_5_all_aug22_m500_test_new_nearest.mat'
tag = 'all_aug22_m500_test_new_nearest' 
fn('/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/rank_k_comparison/rank_k_comparison_stats_m_500_tau_30_eps_5_all_aug22_m500_new_nearest_neighbourconvex.mat')
tag_test_params = 'num_methods=3;m_vals=[500];k_vals=[0.1:0.05:0.3];rho_vals=[0.1:0.1:0.7];phase_type=1;num_trials=50';



load(fn);
eval(tag_test_params);

%% Convert over

    method_labels={};iter=1;
    for new_name = names

        %method_labels{iter}=char(join(new_name{1}));
        method_labels{iter}=char((new_name{1}));
        iter = iter+1;
    end

% Redefine names as needed:    
    %names = {'Neighbour','NuclearNorm'}%,'LRS','BNMF'};%'Convex',
    %names = {'1bit','1bitR','1bitRT'};
    names ={'blank','NIHT','blank'}
    %names = {'TBMC','TBMC-AM','Split-avg','Split-avg-AM','Split-pivot','Split-pivot-AM','Split-bnmf','Split-bnmf-AM','Neighbour','NuclearNorm','1bit','LRS','BNMF'};%'Convex',
mnames=containers.Map(method_labels,names);
    
%% Plots 
%format
%This is a very fixed way of setting the formatting, be careful if you change any of the methods.

if 1
    %This is a very fixed way of setting the formatting, be careful if you change any of the methods.
    
    method_idx=[1:num_methods];
    midx=containers.Map(method_labels,method_idx);
    
    %Now let's set the maps for the different error metrics
    results= {'RESULTS.err_test', 'RESULTS.err_train','RESULTS.comp_time','RESULTS.col_cluster_strength','RESULTS.row_cluster_strength','RESULTS.num_rows_recovered_95','RESULTS.num_rows_recovered_80','RESULTS.num_rows_recovered_75','RESULTS.err_hamming',};
    result_titles = {'Test Error', 'Training Error', 'Computational Time', 'Cluster homogenity','Per cluster error','Proportion of rows recovered to 95\% accuracy','Proportion of rows recovered to 80\% accuracy','Proportion of rows recovered to 75\% accuracy','Overall Error'};
    result_titles_save = {'Test Error', 'Training Error', 'Computational Time', 'Column Recovery','Row Recovery','Proportion of rows recovered to 95 accuracy','Proportion of rows recovered to 80 accuracy','Proportion of rows recovered to 75 accuracy','Overall Error'};
    result_names_y = {'Proportion of entries missclassified','Proportion of entries missclassified','Time, log(s)', 'Strength','Strength','Proportion','Proportion','Proportion','Error'};

    
    mresult_names_y = containers.Map(results,result_names_y);
    mresult_titles = containers.Map(results,result_titles);
    mresult_titles_save = containers.Map(results,result_titles_save);
end


%Script for plotting the phase_transitions for different 
method_idx=[2];%choose the subset you'd like to plot
result_idx=[1:length(results)];
%neighbour
k_plot=[1:length(k_vals)];rho_plot = [1:length(rho_vals)];
%1bitquick
%k_plot = [1:5];rho_plot = [1:4];

top_vals = [ 0.5, 0.5,0,0,0,0,0,0,0.5];
%top_vals =    [ 0.5274, 0.4172,0,0,0,0,0,0,0.4447]; %<- this doesn't
%account for BNMF and NIHT
    for r_idx = [1,2,9]
        result=results{r_idx}        
        close('all');figure('units','normalized','outerposition',[0 0 1 1]);
        plot_val=eval(result);
        pv= plot_val(:,method_idx,k_plot,rho_plot);
        %delete this line; pv= plot_val(:,k_plot,rho_plot);
        %bottom = min(pv(:))
        %top = max(pv(:))
        bottom = 0 ;
        top = top_vals(r_idx);
        top_vals(r_idx) = max(top_vals(r_idx),max(pv(:)));
        [num_trials,~,~,~]=size(plot_val);
        if num_trials>1
            plot_val=squeeze(mean(plot_val));
        else
            plot_val=squeeze(plot_val);
        end
        y=k_vals(k_plot);
        x=rho_vals(rho_plot);
        x_val_label='$1-\rho$';y_val_label='$\phi$';
        for m_idx = method_idx
            close('all')
            method=method_labels{m_idx};
            imagesc(rho_vals(rho_plot),k_vals(k_plot),squeeze(plot_val(midx(method),k_plot,rho_plot)));
            %delete this line!            imagesc(rho_vals,k_vals,squeeze(mean(plot_val)));
            shading interp;
            caxis manual
            caxis([bottom top]);
            set(gca,'YDir','normal','Fontsize',10);
            step_val_x=x(2)-x(1);step_val_y=y(2)-y(1);
            xticks(x);yticks(y);xlim([x(1)-step_val_x/2 x(end)+step_val_x/2]);ylim([y(1)-step_val_y/2 y(end)+step_val_y/2]);
            xlabel(x_val_label,'Interpreter','Latex');
            title_text = string(mresult_titles(result)) +' \\ ' + string(mnames(method))';
            title(sprintf('\\begin{tabular}{c} %s \\end{tabular}',title_text),'interpreter','latex');
            ylabel(y_val_label,'Interpreter','Latex');
            colorbar('TickLabelInterpreter','latex')
            set(gca,'TickLabelInterpreter', 'latex','Fontsize',13)
            saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/ch5plots/tmp/rank_k_comparison_phase_rhok_%s_%s_%s_%s.pdf',result,tag,mresult_titles_save(result),mnames(method)),'pdf')
        end


    end
    
    
    for r_idx=[4,5,6,7,8]     
        result=results{r_idx};
        close('all');figure('units','normalized','outerposition',[0 0 1 1]);
        plot_val=eval(result);
        [num_trials,~,~,~]=size(plot_val);
        if num_trials>1
            plot_val=squeeze(mean(plot_val));
        else
            plot_val=squeeze(plot_val);
        end
        y=k_vals(k_plot);
        x=rho_vals(rho_plot);
        x_val_label='$1-\rho$';y_val_label='$\phi$';
        if  strcmp(result,'comp_time_results')
            %don't try to plot comp time
        else
        for m_idx = method_idx
            close('all')
            method=method_labels{m_idx};
            imagesc(rho_vals(rho_plot),k_vals(k_plot),squeeze(plot_val(midx(method),k_plot,rho_plot)));
            %delete this line:%imagesc(rho_vals,k_vals,squeeze(mean(plot_val)));
            set(gca,'YDir','normal','Fontsize',10);
            xticks(x);yticks(y);            
            step_val_x=x(2)-x(1);step_val_y=y(2)-y(1);
            xticks(x);yticks(y);xlim([x(1)-step_val_x/2 x(end)+step_val_x/2]);ylim([y(1)-step_val_y/2 y(end)+step_val_y/2]);
            xlabel(x_val_label,'Interpreter','Latex');
            title_text = string(mresult_titles(result)) +' \\ ' + string(mnames(method))';
            title(sprintf('\\begin{tabular}{c} %s \\end{tabular}',title_text),'interpreter','latex');
            ylabel(y_val_label,'Interpreter','Latex');
            colorbar('TickLabelInterpreter','latex')
            set(gca,'TickLabelInterpreter', 'latex','Fontsize',13)
            %shg
            saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/ch5plots/tmp/rank_k_comparison_phase_rhok_%s_%s_%s_%s.pdf',result,tag,mresult_titles_save(result),mnames(method)),'pdf')
            %saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/ch5plots/tmp/rank_k_comparison_phase_rhok_%s_%s_%s_%s.eps',result,tag,mresult_titles_save(result),mnames(method)),'epsc')
            %saveas(gcf,sprintf('/home/user/Documents/Mel/Ethera/FirstYear/Images/ch5plots/tmp/rank_k_comparison_phase_rhok_%s_%s_%s_%s.svg',result,tag,mresult_titles_save(result),mnames(method)),'svg')
        end
        end
        
         
    end
    
    
    



    
