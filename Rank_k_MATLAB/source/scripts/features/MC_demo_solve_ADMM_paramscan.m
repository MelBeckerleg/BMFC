%MC_DEMO_SOLVE_ADMM     Solve matrix completion on graphs with ADMM
%
%
%
% see also: (/source/features/utils) MC_solve_ADMM, MC_demo_grid_search, split_observed,
% %           sample_sparse, sample_sparse_t, sample_sparse_AtA, 
%            (/utils) vec
%
%
%code author: Vassilis Kalofolias
%date: Nov 2014
%code edits: Mel Beckerleg
%date: August 2019
%%

close('all');clear;clc
%impath='/home/beckerleg/' ;%
%homepath='/home/beckerleg/Ethera/' ;%
impath='/home/user/Documents/Mel/Ethera/FirstYear/Images/' ;%
homepath='/home/user/Documents/Mel/Ethera/' ;%
addpath([homepath 'BMFC/Rank_k_MATLAB/source/utils'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/ADMM/utils/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/solvers/ADMM/'])
addpath([homepath 'BMFC/Rank_k_MATLAB/source/error/'])
f=figure;g=figure;h1=figure;h2=figure;all_legend=figure;legend1_val=1;legend2_val=1;legend3_val=1;legend4_val=1;legend5_val=1;

for prot_fts={'cathids','none'}%,'funfams','sequences','none'}
    for cmp_fts={'fingerprints','none'}%'fingerprints'
        set_up=1;
        m=200;n=200;
        if set_up 
            load([homepath 'MP2/Data/read_ch/mats/ch_dense_set'],'supported_ch','ch_fingerprints','ch_protein_funfams','ch_protein_sigs'); 
            load([homepath 'MP2/Data/read_ch/mats/ch_protein_cathids.txt'])
            [M,N]=size(supported_ch);
            if 1
                m_idx=randi(M,[m,1]);n_idx=randi(N,[n,1]);
            else
                m_idx=[1:m];n_idx=[1:n];
            end
            
            %Gu=squareform(pdist(ch_protein_funfams(:,[1:3000])'));
            if strcmp(prot_fts,'cathids')
                Gm=ch_protein_cathids(n_idx,:)*ch_protein_cathids(n_idx,:)';
                Gm(find(sum(Gm)==0),find(sum(Gm)==0))=1;
            elseif strcmp(prot_fts,'funfams')
                Gm=ch_protein_funfams(n_idx,:)*ch_protein_funfams(n_idx,:)';
                Gm(find(sum(Gm)==0),find(sum(Gm)==0))=1;
            elseif strcmp(prot_fts,'sequences')
                %this is default; 'none' weights to zero later
                Gm=ch_protein_sigs(n_idx,:)*ch_protein_sigs(n_idx,:)';
                Gm(find(sum(Gm)==0),find(sum(Gm)==0))=1;
            end 
                Gm=Gm./max(Gm);

            if strcmp(cmp_fts,'fingerprints')
                %Gu=squareform(pdist(ch_fingerprints([1:m],:)));
                Gu=ch_fingerprints(m_idx,:)*ch_fingerprints(m_idx,:)';
                Gu(find(sum(Gu)==0),find(sum(Gu)==0))=1;
                Gu=Gu./max(Gu);
            else
                %weighted to zero later!
                Gu=diag(ones(m,1));
            end

            Xn=supported_ch(m_idx,n_idx);
            %%%%%%%%%%%%%%%!!!!!!!!!!!!!!TODO: assign other attributes to struct
            %laplacin= D-A
            Gr=struct;Gr.L=diag(sum(Gu))-Gu;Gr.lmax=max(vec(Gu));
            Gc=struct;Gc.L=diag(sum(Gm))-Gm;Gc.lmax=max(vec(Gm));
            %check

            %% Keep 75% for training, the rest for validation
            rho=0.7;
            [y_train, mask_train, y_val, mask_val, y_test, mask_test] = split_observed(Xn, [rho, (1-rho)/2, (1-rho)/2],0);
            y_train=(y_train+1)/2;y_test=(y_test+1)/2;y_val=(y_val+1)/2;
            Xn((Xn==-1))=0;
            params.size_X = size(Xn);



            %if ~isfield(params, 'zero_mean'), params.zero_mean = 1; end     % this should be true for nuclear norm in general!!
            %if ~isfield(params, 'maxit'), params.maxit = 50; end         % how many iterations?
            %if ~isfield(params, 'verbose'), params.verbose = 0; end
            %if ~isfield(params, 'single'), params.single = isa(y_train, 'single'); end


            %% Normalize data to zero mean and keep the linear transformation details
            y_lims_init = [min(y_train), max(y_train)];

            %mean_train = mean(y_train);
            mean_train=0;
            y_train = y_train - mean_train;   
            y_val = y_val - mean_train;
            %y_test = y_test - mean_train;

            y_lims_scaled = [min(y_train), max(y_train)];

            %% PREPARE PROBLEM PARAMS
            % GRAPHS: (normalized)
            prob_params.Lc = (single(full(Gc.L))/Gc.lmax);
            prob_params.Lr = (single(full(Gr.L))/Gr.lmax);

            %prob_params.Gc_lmax = 1;
            %prob_params.Gr_lmax = 1;

            % DATASETS and masks:
            prob_params.size_X = params.size_X;
            prob_params.mask_val = mask_val;
            prob_params.mask_test = mask_test;
            prob_params.mask_train=mask_train; %added Mel Beckerleg 13/11/2019
            prob_params.A_op = @(x) sample_sparse(x, mask_train);
            prob_params.At_op = @(x) sample_sparse_t(x, mask_train);
            prob_params.AtA_op = @(x) sample_sparse_AtA(x, mask_train);


            %% SOLVER PARAMETERS 
            solver_params.maxit = 200;
            solver_params.verbose = 3;

            solver_params.tol_abs = 2e-6;
            solver_params.tol_rel = 1e-5;

            % need the scaling used for preprocessing to calculate error correctly
            solver_params.y_lims_init = y_lims_init;
            solver_params.y_lims_scaled = y_lims_scaled;

            % MOST IMPORTANT: use verbose = 1 to set rho accordingly (depends on tolerances)
            %solver_params.rho_ADMM = .005000;
            %solver_params.rho_ADMM = .2 * geomean([max(1e-3,prob_params.gamma_n), geomean([max(1e-3,norm(y_train)), max(1e-3,prob_params.gamma_r), max(1e-3,prob_params.gamma_c)])]);

            % for small matrices use false!
            solver_params.svds = false;

        end


        %% Solve the problem using graphs
        %Parameter sweep
   %     error=zeros(15,15);ii=1;jj=1
     %   for weight_val= linspace(0,1,25)
      %      for rho_val=linspace(0,1,25)
            %c_n=1/(150*100);prob_params.gamma_n = c_n*m/2*rho;
            prob_params.gamma_n=0.1;
            if strcmp(cmp_fts,'none')
                prob_params.gamma_r = 0.0;%0.01; %fingerprints is good at about  1/10, generally works well 
                
            else                
                %c_r=1/(10*rho);prob_params.gamma_r = c_r*rho/2;%0.05;%0.05;%0.01;
                prob_params.gamma_r=0.05;
            end
            
            if strcmp(prot_fts,'none')    
                prob_params.gamma_c = 0.0;%0.01;
            else            
                %c_c=1/(10*rho);prob_params.gamma_c = c_c*rho/2;%0.05;%0.01;%0.01; 
                prob_params.gamma_c=0.05;
            end
            solver_params.rho_ADMM = 0.1;%0.009;%rho_val;

            [X_MC_graphs, stat_MC_graphs] = MC_solve_ADMM(y_train, y_val, y_test, prob_params, solver_params);
        %    temp=stat_MC_graphs.rmse_val_round(~isnan(stat_MC_graphs.rmse_val_round));
       %     error(ii,jj)=temp(end)
         %   ii=ii+1
          %  end
           % ii=1
            %jj=jj+1
        %end
        %figure;
        %plot(error)
        %fn=sprintf('/home/beckerleg/Ethera/FirstYear/Images/cathid_scan.eps')
        %saveas(gcf,fn,'epsc')



        %% Solve without graphs, just low rank information
        %prob_params.gamma_n = 3;
        prob_params.gamma_r = 0;
        prob_params.gamma_c = 0;
        %solver_params.rho_ADMM = .15;

        % Now ADMM is equivalent to forward backward algorithm!
        [X_MC_low_rank, stat_MC_low_rank] = MC_solve_ADMM(y_train, y_val, y_test, prob_params, solver_params);
        
        fn=[homepath sprintf( 'stat_MC_regularised_cmp%s_prot%s_small_gammn%s_gamr%s_gamc%s_rho%s_m%s_n%s.mat',string(cmp_fts),string(prot_fts),...
            string(prob_params.gamma_n),string(prob_params.gamma_r),string(prob_params.gamma_c),string(solver_params.rho_ADMM),string(m),string(n))];
        save(fn,'X_MC_low_rank', 'stat_MC_low_rank','X_MC_graphs', 'stat_MC_graphs')
        %% 


        figure(f);
        plot(stat_MC_low_rank.rmse_val);
        hold on
        plot(stat_MC_graphs.rmse_val,'d');
        lgd1{legend1_val}=sprintf('(Benchmark) %s %s',string(cmp_fts),string(prot_fts))
        legend1_val=legend1_val+1;
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
        plot(stat_MC_low_rank.auc_val);
        hold on
        plot(stat_MC_graphs.auc_val,'d');
        %curr=[prev sprintf('lr_cmp%s_prot%s',string(cmp_fts),string(prot_fts)) sprintf('lrg_cmp%s_prot%s',string(cmp_fts),string(prot_fts))]
        %legend(h,curr)
        lgd2{legend2_val}=sprintf('(Benchmark) %s %s',string(cmp_fts),string(prot_fts));
        legend2_val=legend2_val+1;
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
        plot(stat_MC_low_rank.aucpr_val);
        hold on
        plot(stat_MC_graphs.aucpr_val,'d');
        %curr=[prev sprintf('lr_cmp%s_prot%s',string(cmp_fts),string(prot_fts)) sprintf('lrg_cmp%s_prot%s',string(cmp_fts),string(prot_fts))]
        %legend(h,curr)
        lgd3{legend3_val}=sprintf('(Benchmark) %s %s',string(cmp_fts),string(prot_fts));
        legend3_val=legend3_val+1;
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
        p1=plot(stat_MC_low_rank.rmse_val_round);
        hold on
        p2=plot(stat_MC_graphs.rmse_val_round,'d');
        lgd4{legend4_val}=sprintf('(Benchmark) %s %s',string(cmp_fts),string(prot_fts));
        legend4_val=legend4_val+1;
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
        p1=plot(stat_MC_low_rank.rmse_val_round);
        hold on
        p2=plot(stat_MC_graphs.rmse_val_round,'d');
        lgd5{legend5_val}=sprintf('(Benchmark) %s %s',string(cmp_fts),string(prot_fts));
        legend5_val=legend5_val+1;
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



        %% TO PLOT RECOVERED AND INITIAL MATRICES:
        %figure; imagesc(Xn); title('Ground truth (close to low rank with community structure')
        %figure; imagesc(reshape(prob_params.At_op(y_train), params.size_X)); title('20% observations used for recovery')
        %figure; imagesc(X_MC_graphs); title('Recovered matrix from 20% observations, using graph and low rank information')

        %% To reshape in matrix form use:
        % Y_train =   reshape(prob_params.At_op(y_train), params.size_X);
        % Y_val =     reshape(sample_sparse_t(y_val, mask_val), params.size_X);
        % Y_test =    reshape(sample_sparse_t(y_test, mask_test), params.size_X);
    end
end

