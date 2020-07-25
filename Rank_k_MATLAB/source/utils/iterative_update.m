function [x,y] = iterative_update(x,y,W_omega)
    max_iters=20;
    convergence = 0;
    iter = 0;
    while (convergence ~= 1) && (iter<max_iters) 
        x=(W_omega*y>=0)*1.;
        y=((x'*W_omega)>=0)'*1.;
        iter=iter+1;
    end
end