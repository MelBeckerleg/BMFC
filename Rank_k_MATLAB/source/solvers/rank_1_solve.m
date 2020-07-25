function [x,y]= rank_1_solve(W_omega,method)
[m,n]=size(W_omega);
x=zeros(m,1);y=zeros(n,1);
addpath(genpath('/home/user/Documents/Mel/Ethera/Rank_k'))
if strcmp(method,'lp')
    numzs=length(find(W_omega(:)<0));
    if m>1 
        c = -[-ones(numzs,1); 1/2*sum((W_omega')>0)'; 1/2*sum((W_omega)>0)'];
    else 
            c = -[-ones(numzs,1); 1/2*sum((W_omega')>0)'; 1/2*((W_omega)>0)'];
    end

        
    
    %%%this is not the most memory efficient thing to do, but since we need to go up to rho=1 it'll do...
    %%%%we index z such that zk=zij -> zk+1=zi+1j unless i=m.
    %%%%then A is stacking n copies of eye(m) followed by m copies of n
    %%%%ones: [Z,X,Y]
    %%%%%1     1  1
    %     1     1 1
    %      1   1   1
    %       1   1  1  
    %A = [-eye(n*m) kron(ones(n,1),eye(m)) kron(eye(n),ones(m,1)) ];    
    %A=A([find(W_omega(:)<0)],:);
    %A=A(:,[find(W_omega(:)<0)' (m*n+1:1:n*m+m+n)  ]);
    
    A=sparse(numzs,numzs+n+m);
    [I,J]=find(W_omega<0);
    %check to account for the weirdness that happens when W_omega has only
    %one row...
    if m==1
        x=x';y=y';
    end
    A(speye(numzs,numzs)>0)=-ones(numzs,1);
    for ii=1:numzs
        A(ii,numzs+I(ii))=1;
        A(ii,numzs+m+J(ii))=1;
    end
    %A(sub2ind(size(A),(1:numzs)',x+numzs))=1;
    %A(sub2ind(size(A),(1:numzs)',y+m+numzs))=1;
    
   b = [sparse(ones(numzs,1))];

    % Get default options for mosek
    opt = mskoptimset('');
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
    x = out(numzs+1:numzs+m,1);
    y = out(numzs+m+1:end,1);
    
    
    
end

if strcmp(method,'ip')
    numneg=length(find(W_omega(:)<0));numpos=length(find(W_omega(:)>0));
    c=[ones(numneg,1); -ones(numpos,1); zeros(m,1); zeros(n,1)];
    A = [-eye(n*m) kron(ones(m,1),eye(n)) kron(eye(m),ones(n,1)) ;2*eye(n*m) -kron(ones(m,1),eye(n)) -kron(eye(m),ones(n,1))  ];
    A=[A([find(W_omega(:)<0)],:);A([m*n+find(W_omega(:)>0)],:)];
    A=A(:,[find(W_omega(:)<0)' find(W_omega(:)>0)' (m*n+1:1:n*m+m+n)]);
    b = [ones(numneg,1);zeros(numpos,1)];
    options=optimoptions('intlinprog');%options.CutMaxIterations=1;
    out = intlinprog(c,1:length(c),A,b,[],[],zeros(numneg+numpos+m+n,1),ones(numneg+numpos+m+n,1),options);
    %Z = reshape(out(1:numzs,1),m,n);
    x = out(numneg+numpos+1:numneg+numpos+m,1);
    y = out(numneg+numpos+m+1:end,1);  
end

if strcmp(method,'avg')
    x=1.*(sum(W_omega')>0)';
    y=1.*(sum(W_omega)>0)';  
    [m,n]=size(x);
    if (m==1)
        pos_neg_decision=(sum(W_omega)>0);
        y=1.*((W_omega)>-pos_neg_decision)';
    end
end

   
if strcmp(method,'partition')
    idx=randi(n,1);
    x=(W_omega(:,idx)>0)*1;
    if sum(x)
        sub=W_omega(x>0,:);  %
        [msub,nsub]=size(sub);
        if msub>1
            y=(sum(sub)>0)'*1;
        else 
            if (nsub==n) && (sum(sub>0)>0)
                y=(sub>0)'*1.;
            else
                y=zeros(n,1);
            end        
        end
    else
        x=zeros(m,1);
        y=ones(n,1);
    end
end


if strcmp(method,'nmf')
    A=(W_omega+1)/2;
    A(A==1/2)=NaN;
    [W,H]=wnmfrule(A,1);H=H';
    [W,H]=binary_rescale(W,H);
    x=(W>0.5)*1.;
    y=(H>0.5)*1.;
end




end