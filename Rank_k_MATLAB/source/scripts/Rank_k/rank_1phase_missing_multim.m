%fn_home_im='/home/user/Documents/Mel/Ethera/FirstYear/Images/';fn_home='/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/';
%fn_home_im='/home/beckerleg/Ethera/FirstYear/Images/';fn_home='/home/beckerleg/Ethera/FirstYear/Clusters/SecondYear/PhaseWork/';
addpath(genpath('/home/beckerleg/mosek/9.2/toolbox/'));
fn_home_im='/home/beckerleg/Ethera/FirstYear/Images/';fn_home='/home/beckerleg/Ethera/Results/phase_work';
tag='thesis'
bin_thresh=0.05; %Ideally this would be 0
num_trials=50;
step_val=0.05;tau_vals=step_val:step_val:1;tv=length(tau_vals);
rho_vals=0:0.1:1;rv=length(rho_vals);

%parfor m_idx=[3,4,5]
m=300;n=m;
for epsilon=[0.2,0.25,0.3]  
    objective=zeros(tv,tv,rv);error=zeros(tv,tv,rv);error_bin=zeros(tv,tv,rv);error_bin_thresh=zeros(tv,tv,rv);
    l=1;
   for rho=rho_vals
		i=1;j=1;
		for tau1=tau_vals
		    for tau2=tau_vals
                for rand_trials=1:num_trials

                    x_true=zeros(m,1);y_true=zeros(n,1);x_true(1:ceil(tau1*m))=1;y_true(1:ceil(tau2*n))=1; 
                    X=x_true*y_true';
                    X = X + (1-2*X).*(rand(m,n)<epsilon); 

                    W = 2*X-1;
                    W_omega=( W-(1-2*X).*(rand(m,n)<rho));
                    
                    
                    numzs=length(find(W_omega(:)<0));
                    
                    optimal=[zeros(numzs,1);x_true;y_true];
                    %need to set z=1 wherever (5*x_true*y_true'-W_omega)==6; 
                    %and identify the index of z this corresponds to 
                    fi=[find((5*x_true*y_true'-W_omega)==6)];
                    for flipped_idx=1:length(fi)
                        flipped=fi(flipped_idx);
                        optimal(find(find(W_omega<0)==flipped))=1;
                    end



                    
                   
                    c = sparse(-[-ones(numzs,1); 1/2*sum((W_omega')>0)'; 1/2*sum((W_omega)>0)']);
                    %%%this is not the most memory efficient thing to do, but since we need to go up to rho=1 it'll do...
                    %%%%we index z such that zk=zij -> zk+1=zi+1j unless i=m.
                    %%%%then A is  stacking n copied of eye(m) followed by m copies of n ones. 
                    %%%%%1     1  1
                    %     1     1 1
                    %      1   1   1
                    %       1   1  1  
                    A = [-speye(n*m) kron(sparse(ones(m,1)),speye(n)) kron(speye(m),sparse(ones(n,1)))];
                    %check indexing!
                    A=A([find(W_omega(:)<0)],:);
                    A=A(:,[find(W_omega(:)<0)' (m*n+1:1:n*m+m+n)]);
                    b = [sparse(ones(numzs,1))];
                    
                    % Get default options for mosek
                    opt = mskoptimset('')
                    % using simplex to ensure a vertex solution
                    opt = mskoptimset(opt, 'MSK_IPAR_OPTIMIZER','MSK_OPTIMIZER_FREE_SIMPLEX');
                    % using simplex to ensure a vertex solution
                    opt = mskoptimset(opt, 'DISPLAY','OFF');
                   
                    % Set a MOSEK option, in this case turn basic identification off.
                    opt = mskoptimset(opt,'MSK_IPAR_INTPNT_BASIS','MSK_OFF');
                    % Modify a MOSEK parameter with double value
                    opt = mskoptimset(opt,'MSK_DPAR_INTPNT_TOL_INFEAS',1e-12);
                    
                    %without mosek:
                    %opt = optimoptions('linprog','Algorithm','dual-simplex');
                    out = linprog(c,A,b,[],[],zeros(numzs+m+n,1),ones(numzs+m+n,1),opt);
                    %Z = reshape(out(1:numzs,1),m,n);
                    u = out(numzs+1:numzs+m,1);
                    v = out(numzs+m+1:end,1);

                    %success?
                    error_bin(i,j,l)=error_bin(i,j,l)+1-(((norm(u-x_true)+norm(v-y_true)))>0)*1.;
                    error_bin_thresh(i,j,l)=error_bin_thresh(i,j,l)+1-(((norm(u-x_true)^2+norm(v-y_true)^2)/(m+n))>bin_thresh)*1.;
                    error(i,j,l)=error(i,j,l)+(sum(abs(u-x_true))+sum(abs(v-y_true)))/(n+m);
                    c'*(out-optimal);
                    objective(i,j,l)=objective(i,j,l)+c'*(out-optimal);
                    fn=sprintf('%serror_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
                    save(fn,'error')
                    fn=sprintf('%serrorbin_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
                    save(fn,'error_bin')
                    fn=sprintf('%serrorbinthresh_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
                    save(fn,'error_bin_thresh')
                    fn=sprintf('%serrorobjective_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
                    save(fn,'objective')
                end
                error(i,j,l)=error(i,j,l)/num_trials;
                error_bin(i,j,l)=error_bin(i,j,l)/num_trials;
                error_bin_thresh(i,j,l)=error_bin_thresh(i,j,l)/num_trials;
                j=j+1;
		    end
		    i=i+1; j=1;
        end
        i=1;j=1;l=l+1;
	end

end
%end
%%%%%%%%%%%%Plots




 %% Plots on my computer
if 0
    m=100;
    
    %confirmation 
    %step_val=0.1
    %rho_vals=0:0.1:0.1;
    
    %thesis
    step_val=0.05;
    rho_vals=0:0.1:1;
    
    %all
   
    tau_vals=0.1:step_val:1;tau_vals=tau_vals(1:end-1);
    
    
    fn_home_im='/home/user/Documents/Mel/Ethera/FirstYear/Images/tmp'
    
    %confirmation: fn_home='/home/user/Documents/Mel/Ethera/Results/'%folder temp for ntrials=3
    %thesis: 
    fn_home='/home/user/Documents/Mel/Ethera/BMFC/Rank_k_MATLAB/Results/phase_work';
    set(0,'defaulttextinterpreter','latex')

    for epsilon=[0.2,0.25,0.3]
        
        %confirmation:
        
%         fn=sprintf('%serror_epsilon%s',fn_home,string(100*epsilon));
%         load(fn,'error')
%         fn=sprintf('%serrorbin_epsilon%s',fn_home,string(100*epsilon));
%         load(fn,'error_bin')
%         fn=sprintf('%serrorbinthresh_epsilon%s',fn_home,string(100*epsilon));
%         load(fn,'error_bin_thresh')
%         fn=sprintf('%serrorobjective_epsilon%s',fn_home,string(100*epsilon));
%         load(fn,'objective')
        
        %thesis:
        fn=sprintf('%serror_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
        load(fn,'error')
        fn=sprintf('%serrorbin_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
        load(fn,'error_bin')
        fn=sprintf('%serrorbinthresh_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
        load(fn,'error_bin_thresh')
        fn=sprintf('%serrorobjective_epsilon%s%s%s',fn_home,string(100*epsilon),string(m),tag);
        load(fn,'objective')


        for error_type={'recovered','hamming','threshrecovered','objective'}
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
        
            
        for rho_val=1:3
            rho=rho_vals(rho_val);
            figure();

            imagesc(plot_var(1:end-1,1:end-1,rho_val))
            set(gca,'YDir','normal','Fontsize',20);xlabel('$\tau_1$','Interpreter','Latex');ylabel('$\tau_2$','Interpreter','Latex')
            
            xticklabels=[0.1,0.3,0.5,0.7,0.9];
            xticks=linspace(1,size(plot_var(1:end-1,1:end-1,1),2),numel(xticklabels));
            
            set(gca,'Xtick', xticks-step_val,'XTickLabel',xticklabels,'Ytick',xticks,'YTickLabel',xticklabels,'TickLabelInterpreter', 'latex')
            %check calibration: hold on; plot(1,1,'*');plot(9,9,'*')
            fun1=@(a,e) e*(a-1)./((-1+4*a)*e-a);
            %rescale
            plotscale=@(truth) truth/step_val;%(truth-0.1/2)/step_val;%(truth-tau_vals(1))*1/step_val+1;
            y=zeros(length(tau_vals),1);i=0;y1=y;y2=y;y3=y;
            y1(1)=plotscale(fun1(0,epsilon));i=2;
            for a=tau_vals                
                y1(i)=plotscale(fun1(a,epsilon));
                i=i+1;
            end
            hold on; plot([plotscale(0),plotscale(tau_vals)],y1,'color','k','linewidth',2,'linestyle','--');
            val=plotscale(epsilon/(2*(1-epsilon)));
            line([0,numel(tau_vals)],[val,val],'color','k','linewidth',2,'linestyle',':');line([val,val],[0,numel(tau_vals)],'color','k','linewidth',2,'linestyle',':');
            xlim([plotscale(step_val/2),numel(tau_vals)]);ylim([plotscale(step_val/2),numel(tau_vals)]);colorbar
            set(gcf, 'PaperPositionMode', 'auto');colorbar('TickLabelInterpreter','latex')
            fn=sprintf('%stau1tau2recoveryphaseepsilon%s%sjune%s',fn_home_im,string(100*epsilon),string(error_type),string(rho_val));
            saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close

        end
            
 
            
         

        %Plot how recovery changes with rho for different Tau. 
            markers={'^','o','+','>','<','*','s','d','p',''};
            set(0,'defaulttextInterpreter','latex');set(legend,'Interpreter','latex')
            figure()
            for i=1:2:length(tau_vals)
                hold on 
            legend('AutoUpdate','on')
            plot(rho_vals,reshape(plot_var(i,i,:),1,[]),sprintf('-%s',markers{ceil(i/2)}),'DisplayName',sprintf('Tau=%0g\n',round(tau_vals(i),2,'significant')),'LineWidth',2,'MarkerSize',10)
            legend('AutoUpdate','on')
            end
            hold off

            xlabel('$\rho$')
            ylabel(sprintf('Error:%s',string(error_type)))
            fn=sprintf('%stau1tau2recoveryphasemissingepsilon%s%sjune',fn_home_im,string(100*epsilon),string(error_type));
            saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close
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
end