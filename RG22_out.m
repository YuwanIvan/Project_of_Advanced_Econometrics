function [h,z,u]=RG22_out(theta,ret,x)
    
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

    
    % initial values
    h00=var(ret);
    h0=var(ret);
    x00=mean(x);
    x0=mean(x);

    
    h=zeros(T,1);
    
    h(1)=exp(w+b1*log(h0)+b2*log(h00)+r1*log(x0)+r2*log(x00));
    h(2)=exp(w+b1*log(h(1))+b2*log(h0)+r1*log(x(1))+r2*log(x0));

    for t=3:T
        h(t)=exp(w+b1*log(h(t-1))+b2*log(h(t-2))+r1*log(x(t-1))+...
            r2*log(x(t-2)));
    end

    z=ret./sqrt(h);
    u=log(x)-ks-ph*log(h)-tau1*z-tau2*(z.^2-1);

    
end