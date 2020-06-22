m = 100;
n = 100;
num_trials=50;
step_val=0.01;tau_vals=0.1:step_val:1;tv=length(tau_vals);
for epsilon=[0.2,0.25,0.3]
   
error=zeros(length(tau_vals),length(tau_vals));error_bin=zeros(length(tau_vals),length(tau_vals));
i=1;j=1;
for tau1=tau_vals
    for tau2=tau_vals
for rand_trials=1:num_trials

%tau1 = 0.25;
%tau2 = 0.25;

%X = [ones(ceil(tau1*m),ceil(tau2*n)) zeros(ceil(tau1*m),floor((1-tau2)*n)); zeros(floor((1-tau1)*m),n)];
x_true=zeros(m,1);y_true=zeros(n,1);x_true(1:ceil(tau1*m))=1;y_true(1:ceil(tau2*n))=1;
X=x_true*y_true';
X = X + (1-2*X).*(rand(m,n)<epsilon); 

W = 1 - 2*X;

c = [W(:); zeros(m+n,1)];
A = [-eye(m*n) kron(eye(n),ones(m,1)) kron(ones(n,1),eye(m)); 2*eye(m*n) -kron(eye(n),ones(m,1)) -kron(ones(n,1),eye(m))];
b = [ones(m*n,1); zeros(m*n,1)];

options = optimoptions('linprog','Algorithm','dual-simplex');
out = linprog(c,A,b,[],[],zeros(m*n+m+n,1),ones(m*n+m+n,1),options);

Z = reshape(out(1:m*n,1),m,n);
v = out(m*n+1:m*n+m,1);
u = out(m*n+m+1:end,1);
%plot
%imagesc(floor(Z));

%success?
tol=0.1;
error_bin(i,j)=error_bin(i,j)+1-((norm(u-x_true)+norm(v-y_true))>0)*1.;
error(i,j)=error(i,j)+(sum(abs(u-x_true))+sum(abs(v-y_true)))/(n+m);
fn=sprintf('/home/beckerleg/Ethera/FirstYear/Clusters/SecondYear/EvilPlans/error_epsilon%s',string(epsilon));
save(fn,'error')
save('/home/beckerleg/Ethera/FirstYear/Clusters/SecondYear/EvilPlans/error_bin','error_bin')
end
error(i,j)=error(i,j)/num_trials;
error_bin(i,j)=error_bin(i,j)/num_trials;
j=j+1;
    end
    i=i+1;
    j=1;
end

for error_type={'recovered','hamming'}
	if strcmp(error_type,'recovered')
		imagesc(error_bin);

	elseif strcmp(error_type,'hamming')
		imagesc(error);
end
	set(gca,'YDir','normal');xlabel('$\tau_1$','Interpreter','Latex');ylabel('$\tau_2$','Interpreter','Latex')
	xticklabels=tau_vals([1,ceil(length(tau_vals)/4),ceil(length(tau_vals)/2),ceil(3*length(tau_vals)/4),length(tau_vals)]) ;
	xticks=linspace(1,size(error,2),numel(xticklabels));
	set(gca,'Xtick', xticks,'XTickLabel',xticklabels,'Ytick',xticks,'YTickLabel',xticklabels)
	fun1=@(a,e) e*(a-1)./((-1+4*a)*e-a);plotscale=@(truth) (truth-tau_vals(1))*1/step_val+1;
	y=zeros(length(tau_vals),1);i=0;y1=y;y2=y;y3=y;
	for a=tau_vals
	    i=i+1;
	    y1(i)=plotscale(fun1(a,epsilon));
	end
	hold on; plot(y1,'color','k');
	val=plotscale(epsilon/(2*(1-epsilon)));
	line([0,numel(tau_vals)],[val,val]);line([val,val],[0,numel(tau_vals)]);
	xlim([0,numel(tau_vals)]);ylim([0,numel(tau_vals)]);
colorbar


	fn=sprintf('/home/beckerleg/Ethera/FirstYear/Clusters/SecondYear/EvilPlans/tau1tau2recoveryphaseepsilon%s%s',string(epsilon),string(error_type))
	saveas(gca,fn,'epsc');saveas(gca,fn,'pdf');close
	end


end
