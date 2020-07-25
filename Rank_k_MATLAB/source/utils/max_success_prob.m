function [min_prob] = max_success_prob(n,a,rho,k)
t1=(1-a)./(1-a.^k);

min_prob=1;min_i=1;
for ii=1:k-1
    delta=(2-3*a.^2+a.^(2*(k-ii)))./(2-a.^2-a.^(2*(k-ii)));
    %(1-a.^2-a.^(2*(k-ii)))./(3-a.^2-a.^(2*(k-ii)));
    mu1=t1.^2.*n.^2.*rho*a^(2*ii)    
    mu2=1/2*t1.^2.*n.^2.*rho.*a^(2*ii+2).*((1-a.^(2*(k-ii-1)))./(1-a.^2))
    mu2tilde=1/2*t1.^2.*n.^2.*rho.*a^(2*ii+4).*((1-a.^(2*k-2*ii-2))./(1-a.^2))
    probA=1-exp(-delta.^2.*mu1./(2))-exp(-delta.^2.*mu1./(2+delta))
    probB=1-exp(-delta.^2.*mu2./(2+delta))
    probBgivenA=1-exp(-delta^2*mu2tilde/(2+delta))
    prob=probA+probB
    if prob<min_prob
       min_prob=prob;
       min_i=ii;
    end
end

k-min_i


end