function [min_prob] = min_success_prob_independent(n,a,rho,k)
t1=(1-a)./(1-a.^k);

min_prob=1;min_i=1;
for ii=1:k-2
    delta=(1-a.^2-a.^(2*(k-ii)))./(3-a.^2-a.^(2*(k-ii)));
    mu1=1/2*t1.^2.*n.^2.*rho.*a^(2*ii)
    mu2=1/2*t1.^2.*n.^2.*rho./2*a^(2*ii+2).*((1-a.^(k-ii))./(1-a.^2))
    probA=1-exp(-delta.^2.*mu1./(2+delta))
    probB=1-exp(-delta.^2.*mu2./(2))
    prob=probA+probB-1
    if prob<min_prob
       min_prob=prob;
       min_i=ii;
    end
end

k-min_i

end