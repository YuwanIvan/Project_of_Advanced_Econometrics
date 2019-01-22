function logLp=RG22s_logLp(theta,ret,x)
    % This is the partital log-likelihood function of the full
    % model, that is, RealGARCH(2,2) with lagged returns added.
    % The output is a scalar.
    
    T=length(ret);
    
    w=theta(1);
    b1=theta(2);
    b2=theta(3);
    r1=theta(4);
    r2=theta(5);
    ks=theta(6);
    ph=theta(7);
    tau1=theta(8);
    tau2=theta(9);
    sig2=theta(10);
    a=theta(11);
    
    % initial values
    h00=var(ret);
    h0=var(ret);
    x00=mean(x);
    x0=mean(x);
    ret0=mean(ret);
    
    eps=10^(-20);
    
    h=zeros(T,1);
    h(1)=exp(w+b1*log(h0)+b2*log(h00)+r1*log(x0)+r2*log(x00)+a*log(max(ret0^2,eps)));
    h(2)=exp(w+b1*log(h(1))+b2*log(h0)+r1*log(x(1))+r2*log(x0)+a*log(max(ret(1)^2,eps)));
    for t=3:T
        h(t)=exp(w+b1*log(h(t-1))+b2*log(h(t-2))+r1*log(x(t-1))+...
            r2*log(x(t-2))+a*log(max(ret(t-1)^2,eps)));
    end
    
    logLp=sum(-0.5*log(2*pi)-0.5*log(h)-(ret.^2)./(2*h));
    
end