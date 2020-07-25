%Plotting the values after a test.

%% Set the following variables
%string(prob_params.gamma_n),string(prob_params.gamma_r),string(prob_params.gamma_c),string(solver_params.rho_ADMM),string(m),string(n));

%% mapping

%{'cathids','none','funfams','sequences','none'}
%{'fingerprints','none'}

f=figure;g=figure;h1=figure;h2=figure;all_legend=figure;legend1_val=1;legend2_val=1;legend3_val=1;legend4_val=1;legend5_val=1;

for prot_fts={'cathids','none','funfams','sequences','none'}
    for cmp_fts={'fingerprints','none'}%'fingerprints'
        
        %% load in the test
        
        fn=[homepath sprintf( 'stat_MC_regularised_cmp%s_prot%s_small_gammn%s_gamr%s_gamc%s_rho%s_m%s_n%s.mat',string(prob_params.col_features),string(prob_params.row_features),...
        string(prob_params.gamma_n),string(prob_params.gamma_r),string(prob_params.gamma_c),string(solver_params.rho_ADMM),string(m),string(n))];
        load(fn,'X_MC_graphs', 'stat_MC_graphs')

        
        figure(f);
        hold on
        plot(stat_MC_graphs.rmse_val,'d');
        lgd1{legend1_val}=sprintf('%s %s',string(cmp_fts),string(prot_fts))
        legend1_val=legend1_val+1;
        %leg1=legend(lgd1,'interpreter','latex','Location','southoutside');
        title('Validation rmse of different models','Interpreter','latex');
        xlabel('Iteration','Interpreter','latex');
        ylabel('RMSE','Interpreter','latex');
        ax=gca;ax.TickLabelInterpreter='latex';set(gca,'FontSize',15);
        fn=[impath sprintf( 'regularisation_cmp%s_prot%s_small.eps',string(cmp_fts),string(prot_fts))];
        saveas(gcf,fn,'epsc');
        
        

        figure(h1);
        hold on
        plot(stat_MC_graphs.auc_val,'d');
        %curr=[prev sprintf('lr_cmp%s_prot%s',string(cmp_fts),string(prot_fts)) sprintf('lrg_cmp%s_prot%s',string(cmp_fts),string(prot_fts))]
        %legend(h,curr)
        lgd2{legend2_val}=sprintf('%s %s',string(cmp_fts),string(prot_fts));
        legend2_val=legend2_val+1;
        %legend(lgd2,'interpreter','latex','Location','southoutside')
        title('Validation auc of different models','Interpreter','Latex')
        xlabel('iteration','Interpreter','Latex')
        ylabel('AUC','Interpreter','Latex')
        ax=gca;ax.TickLabelInterpreter='latex';set(gca,'FontSize',15);
        fn=[impath sprintf('regularisation_cmp%s_prot%s_small_auc.eps',string(cmp_fts),string(prot_fts))]
        saveas(gcf,fn,'epsc')
        
        figure(h2);
        hold on
        plot(stat_MC_graphs.aucpr_val,'d');
        %curr=[prev sprintf('lr_cmp%s_prot%s',string(cmp_fts),string(prot_fts)) sprintf('lrg_cmp%s_prot%s',string(cmp_fts),string(prot_fts))]
        %legend(h,curr)
        lgd3{legend3_val}=sprintf('%s %s',string(cmp_fts),string(prot_fts));
        legend3_val=legend3_val+1;
        %legend(lgd3,'interpreter','latex','Location','southoutside')
        title('Validation auc of different models','Interpreter','Latex')
        xlabel('iteration','Interpreter','Latex')
        ylabel('AUC-PR','Interpreter','Latex')
        ax=gca;ax.TickLabelInterpreter='latex';set(gca,'FontSize',15);
        fn=[impath sprintf('regularisation_cmp%s_prot%s_small_aucpr.eps',string(cmp_fts),string(prot_fts))]
        saveas(gcf,fn,'epsc')


        figure(g);
        hold on
        p2=plot(stat_MC_graphs.rmse_val_round,'d');
        lgd4{legend4_val}=sprintf('%s %s',string(cmp_fts),string(prot_fts));
        legend4_val=legend4_val+1;
        %legend(lgd3,'interpreter','latex','Location','southoutside')
        title('Rounded validation error of different models','Interpreter','Latex');
        xlabel('iteration','Interpreter','Latex');
        ylabel('RMSE','Interpreter','Latex');
        ax=gca;ax.TickLabelInterpreter='latex';set(gca,'FontSize',15);
        fn=[impath sprintf('rounded_regularisation_cmp%s_prot%s_small.eps',string(cmp_fts),string(prot_fts))]
        saveas(gcf,fn,'epsc')

        
        figure(all_legend);
        p2=plot(stat_MC_graphs.rmse_val_round,'d');
        lgd5{legend5_val}=sprintf('%s %s',string(cmp_fts),string(prot_fts));
        legend5_val=legend5_val+1;
        %legend(lgd3,'interpreter','latex','Location','southoutside')
        title('Rounded validation error of different models','Interpreter','Latex');
        xlabel('iteration','Interpreter','Latex');
        ylabel('RMSE','Interpreter','Latex');
        ax=gca;ax.TickLabelInterpreter='latex';set(gca,'FontSize',15);
        %set(p1,'visible','off');set(p2,'visible','off');set(ax,'visible','off');
        legend(lgd5,'interpreter','latex','Location','southoutside')
        fn=[impath sprintf('LEGEND_rounded_regularisation_cmp%s_prot%s_small.eps',string(cmp_fts),string(prot_fts))]
        saveas(gcf,fn,'epsc')
    end
end