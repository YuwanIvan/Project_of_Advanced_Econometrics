function [h,z,u]=RG11_out(theta,ret,x)
    
    T=length(ret);
    
    w=theta(1);
    b1=theta(2);

    r1=theta(3);

    ks=theta(4);
    ph=theta(5);
    tau1=theta(6);
    tau2=theta(7);
    sig2=theta(8);

    
    % initial values

    h0=var(ret);

    x0=mean(x);

    
    h=zeros(T,1);
    
    h(1)=exp(w+b1*log(h0)+r1*log(x0));

    for t=2:T
        h(t)=exp(w+b1*log(h(t-1))+r1*log(x(t-1)));
    end

    z=ret./sqrt(h);
    u=log(x)-ks-ph*log(h)-tau1*z-tau2*(z.^2-1);

    
end