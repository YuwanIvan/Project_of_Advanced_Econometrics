function logf=G11_logf(theta,ret)
    % This is the observation-specific log-likelihood function of 
    % GARCH(1,1) with lagged returns added.
    % The output is a T-dimensional vector.
    
    T=length(ret);
    
    w=theta(1);
    b1=theta(2);
    a=theta(3);
    
    % initial values
    h0=var(ret);
    ret0=mean(ret);
    
    eps=10^(-20);
    
    h=zeros(T,1);
    
    h(1)=exp(w+b1*log(h0)+a*log(max(ret0^2,eps)));

    for t=2:T
        h(t)=exp(w+b1*log(h(t-1))+a*log(max(ret(t-1)^2,eps)));
    end

    z=ret./sqrt(h);

    
    logf=-0.5*log(2*pi)-0.5*log(h)-(ret.^2)./(2*h);
    
end