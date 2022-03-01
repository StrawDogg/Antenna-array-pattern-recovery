function s=levyFlight()
    beta = 1.5;
    sigma_u = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    sigma_v = 1;
    u = normrnd(0,sigma_u);
    v = normrnd(0,sigma_v);
    s = u/(abs(v))^(1/beta);
end 