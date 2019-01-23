function logf=RG21_logf(theta,ret,x)
    % This is the observation-specific log-likelihood function of
    % RealGARCH(2,1) model.
    % The output is a T-dimensional vector.
    
    T=length(ret);
    
    w=theta(1);
    b1=theta(2);
    b2=theta(3);
    r1=theta(4);

    ks=theta(5);
    ph=theta(6);
    tau1=theta(7);
    tau2=theta(8);
    sig2=theta(9);

    
    % initial values
    h00=var(ret);
    h0=var(ret);

    x0=mean(x);

    
    h=zeros(T,1);
    
    h(1)=exp(w+b1*log(h0)+b2*log(h00)+r1*log(x0));
    h(2)=exp(w+b1*log(h(1))+b2*log(h0)+r1*log(x(1)));

    for t=3:T
        h(t)=exp(w+b1*log(h(t-1))+b2*log(h(t-2))+r1*log(x(t-1)));
    end

    z=ret./sqrt(h);
    u=log(x)-ks-ph*log(h)-tau1*z-tau2*(z.^2-1);
    
    logf=-0.5*log(2*pi)-0.5*log(h)-(ret.^2)./(2*h)...
        -0.5*log(2*pi*sig2)-0.5*(u.^2)./sig2;
    
end
