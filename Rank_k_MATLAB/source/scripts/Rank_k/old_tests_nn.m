%Define a function to perform nearest neighbours classification on
%databases
clear;clc;
addpath('/home/user/Matlab')
addpath('/home/user/Documents/Mel/Ethera/SecondYear/ChernoffBounds/')
%Set-up
%then change mvals to rho_vals
xval='m';
tau=0.3;epsilon=0.2;
num_trials=100;
if(strcmp(xval,'rho'))
    mvalues=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
    m=50;n=50;
    r=3;%floor(sqrt(m));
else
    if(strcmp(xval,'m'))
    mvalues=[100,200,300,400,500];
    rho=0.3;
    else
        if(strcmp(xval,'r'))
            m=100;n=100;
            rho=0.3;
            %mvalues=[1:ceil(sqrt(m)/10):sqrt(m)];
            mvalues=[1:1:10];
        end
    end
end

%
mvals=length(mvalues);

error_nn=zeros(num_trials,mvals);xerror_nn=zeros(num_trials,mvals);yerror_nn=zeros(num_trials,mvals);time_nn=zeros(num_trials,mvals);
error_c=zeros(num_trials,mvals);xerror_c=zeros(num_trials,mvals);yerror_c=zeros(num_trials,mvals);time_c=zeros(num_trials,mvals);
error_cr=zeros(num_trials,mvals);xerror_cr=zeros(num_trials,mvals);yerror_cr=zeros(num_trials,mvals);time_cr=zeros(num_trials,mvals);
error_mc=zeros(num_trials,mvals);

for mval=1:mvals
    if(strcmp(xval,'rho'))
        rho=mvalues(mval);
    else
        if(strcmp(xval,'m'))
        m=mvalues(mval);n=m;
        r=20;
        
        else
            if(strcmp(xval,'r'))
            r=mvalues(mval);
            end
        end   
    end
for trial=1:num_trials

clusters=randperm(m);
X=zeros(m,r);Y=zeros(n,r);
for rank=1:r
    X(clusters((rank-1)*floor(m/r)+1:floor(m/r)*rank), rank)=1;
    Y(:,rank)=(rand(n,1)>1-tau)*1;
end
A=X*Y';
noise_mask=(rand(m,n)<epsilon)*1;
A_noise=A-noise_mask.*(2*A-1);
mask=(rand(m,n)>rho)*1;

tic
[X_approx,Y_approx]=NN(A_noise,mask,r);
time_nn(trial,mval)=toc;
error_nn(trial,mval)=(norm(A-X_approx*Y_approx','fro')/sqrt(m*n))^2;
xerror_nn(trial,mval)=cluster_error(X,X_approx);
yerror_nn(trial,mval)=cluster_error(Y,Y_approx);

tic
[X_approx,Y_approx]=Convex(A_noise,mask,r,0,1);
time_c(trial,mval)=toc;
error_c(trial,mval)=(norm(A-X_approx*Y_approx','fro')/sqrt(m*n))^2;
xerror_c(trial,mval)=cluster_error(X,X_approx);
yerror_c(trial,mval)=cluster_error(Y,Y_approx);

tic
[X_approx,Y_approx]=Convex(A_noise,mask,r,1,1);
time_cr(trial,mval)=toc;
error_cr(trial,mval)=(norm(A-X_approx*Y_approx','fro')/sqrt(m*n))^2;
xerror_cr(trial,mval)=cluster_error(X,X_approx);
yerror_cr(trial,mval)=cluster_error(Y,Y_approx);

%benchmark minimal error
[CompletedMat, ier] = MatrixCompletion(A_noise.*mask, mask,100, 'nuclear', 10, 1e-8, 0);
error_mc(trial,mval)=norm((A-CompletedMat>0.5)*1.);
end
end

tag=sprintf('DEClarge%s%d',xval, r);
%tag=['low_noise',xval];

close all
figure()
set(gca,'Fontsize',25,'TickLabelInterpreter','latex')
ylabel('Error','Interpreter','latex')
xlabel(sprintf('$%s$',xval),'Interpreter','latex')
hold on
errorbar(mvalues,mean(error_nn),std(error_nn),'LineWidth',2,'Marker','*','MarkerSize',10)
errorbar(mvalues,mean(error_c),std(error_c),'LineWidth',2,'Marker','d','MarkerSize',10);
errorbar(mvalues,mean(error_cr),std(error_cr),'LineWidth',2,'Marker','o','MarkerSize',10)
%errorbar(mvalues,mean(error_mc),std(error_mc))
l=legend('NN','Convex','ConvexR','location','NorthWest');
set(l,'Interpreter','latex');
saveas(gcf,sprintf(['/home/user/Documents/Mel/Ethera/FirstYear/Images/rho_total',tag]),'epsc')
figure()
set(gca,'Fontsize',25,'TickLabelInterpreter','latex')
ylabel('Row recovery score','Interpreter','latex')
xlabel(sprintf('$%s$',xval),'Interpreter','latex')
hold on
errorbar(mvalues,mean(xerror_nn),std(xerror_nn),'LineWidth',2,'Marker','*','MarkerSize',10)
errorbar(mvalues,mean(xerror_c),std(xerror_c),'LineWidth',2,'Marker','d','MarkerSize',10);
errorbar(mvalues,mean(xerror_cr),std(xerror_cr),'LineWidth',2)
%l=legend('NN','Convex','ConvexR','location','SouthEast');
set(l,'Interpreter','latex');
saveas(gcf,sprintf(['/home/user/Documents/Mel/Ethera/FirstYear/Images/rho_row',tag]),'epsc')
figure()
set(gca,'Fontsize',25,'TickLabelInterpreter','latex')
ylabel('Column recovery score','Interpreter','latex')
xlabel(sprintf('$%s$',xval),'Interpreter','latex')
hold on
errorbar(mvalues,mean(yerror_nn),std(yerror_nn),'LineWidth',2,'Marker','*','MarkerSize',10)
errorbar(mvalues,mean(yerror_c),std(yerror_c),'LineWidth',2,'Marker','d','MarkerSize',10);
errorbar(mvalues,mean(yerror_cr),std(yerror_cr),'LineWidth',2,'Marker','o','MarkerSize',10);
%l=legend('NN','Convex','ConvexR','location','SouthEast');
set(l,'Interpreter','latex');
saveas(gcf,sprintf(['/home/user/Documents/Mel/Ethera/FirstYear/Images/rho_col',tag]),'epsc')

figure()
set(gca,'Fontsize',25,'TickLabelInterpreter','latex')
ylabel('Time, s','Interpreter','latex')
xlabel(sprintf('$%s$',xval),'Interpreter','latex')
hold on 
errorbar(mvalues,mean(time_nn), std(time_nn),'LineWidth',2,'Marker','*','MarkerSize',10)
errorbar(mvalues,mean(time_c), std(time_nn),'LineWidth',2,'Marker','*','MarkerSize',10)
errorbar(mvalues,mean(time_cr), std(time_nn),'LineWidth',2,'Marker','*','MarkerSize',10)
l=legend('NN','Convex','ConvexR','location','SouthEast');
set(l,'Interpreter','latex');
saveas(gcf,sprintf(['/home/user/Documents/Mel/Ethera/FirstYear/Images/time',tag]),'epsc')

% 
% close all
% figure()
% title('Error')
% hold on
% plot(error_nn)
% plot(error_c)
% plot(error_cr)
% plot(error_MC)
% legend('NN','Convex','ConvexR','MC')
% figure()
% title('Row recovery error')
% hold on
% plot(xerror_nn)
% plot(xerror_c)
% plot(xerror_cr)
% legend('NN','Convex','ConvexR')
% figure()
% title('Column recovery error')
% hold on
% plot(yerror_nn)
% plot(yerror_c)
% plot(yerror_cr)
% legend('NN','Convex','ConvexR')



% function [val] = cluster_error(X,Xapprox,r) 
% val=sum(sum(abs(X)))+sum(sum(abs(Xapprox)));
% permutations=perms(1:r);
% for prank=1:r
%     nval=sum(sum(abs(Xapprox(permutations(prank),:)-X)))/r;
%     if nval<val
%        val=nval;
%     end
% end
% 
% end

function [val] = cluster_error_small_r(X,Xapprox,r) 
val=0;
for rval=1:r
val=val+mean(pdist(X(Xapprox(:,rval)>0,:)));
end
val=1/r*val;
end


function [val] = cluster_error(X,Xapprox) 
[m,r]=size(X);
base_val1=max(max(X'*X));
base_val2=max(max((1-X)'*(1-X)));
val1=1/(base_val1*r)*min(sum(abs(max(X'*Xapprox))),sum(abs(max((X'*Xapprox))')));
val2=1/(base_val2*r)*min(sum(abs(max((1-X)'*(1-Xapprox)))),sum(abs(max(((1-X)'*(1-Xapprox)))')));
val=min(val1,val2);
end
    
%evaluation of the error

function D2 = nanhamdist(XI,XJ)  
%NANHAMDIST Hamming distance ignoring coordinates with NaNs
[m,p] = size(XJ);
nesum = zeros(m,1);
pstar = zeros(m,1);
for q = 1:p
    notnan = ~(isnan(XI(q)) | isnan(XJ(:,q)));
    nesum = nesum + ((XI(q) ~= XJ(:,q)) & notnan);
    pstar = pstar + notnan;
end
D2 = nesum./pstar; 
end

function [X,Y]= NN(A,mask,r)
observed_hamming=0;
A=A.*mask;
[m,n]=size(A);
exact_rho=sum(sum(mask))/m/n;
[m,n]=size(A);
if observed_hamming
    A_copy=A;
    A_copy(mask)=Nan;
    H=squareform(pdist(A_copy,@nanhamdist));
else
    H=squareform(pdist(A,'hamming'));
end
X=zeros(m,r);
Y=zeros(n,r);
hvals=ones(m,1);
    for j=1:r
        Hsub=H(find(hvals),:);
        [val,idx]=sort(Hsub(j,:));
        idx=idx(1:m/r);
        C=idx;
        for i=1:length(C)
            iidx=C(i);
            X(iidx,j)=1;
            hvals(iidx)=0;
        end
        mu=(1/(exact_rho*length(idx))*sum(A(idx,:))>1/2)*1;
        Y(:,j)=mu;
    end
end

%Define a function to do the convex minimisation

function [X,Y]=Convex(A,B,r,recovered,binarised)
[m,n]=size(A);
exact_rho=sum(sum(B))/m/n;
X=zeros(m,r);
Y=zeros(n,r);
[CompletedMat, ier] = MatrixCompletion(A.*B, B,100, 'nuclear', 10, 1e-8, 0);
if binarised
    CompletedMat=(CompletedMat>0.5)*1;
end
H=squareform(pdist(CompletedMat));
labels=zeros(m,1);
label=1;
%find X

for i=1:m
    if labels(i)==0
        [val,idx]=sort(H(i,:));
        idx=idx(1:ceil(m/r/2));
        if sum(labels(idx))
            existing_label=nonzeros(labels(idx));
            labels(idx)=existing_label(1);
        else
            labels(idx)=label;
            label=label+1;
        end
    end
end
unique_labels=unique(labels);
recovered_rank=min(length(unique_labels),r);
for j=1:recovered_rank
    X(find(labels==unique_labels(j)),j)=1;
end

%find Y

if recovered
    mat=CompletedMat;
else
    mat=A.*B;    
end
for j=1:recovered_rank
    mu=(1/(exact_rho*(length(find(labels==unique_labels(j)))))*sum(mat(find(labels==unique_labels(j)),:)))'>1/2;
    Y(:,j)=mu;
end
    
end
    

