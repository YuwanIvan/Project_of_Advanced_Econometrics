function [h,z,u]=RG11lin_out(theta,ret,x)
    
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
    
    h(1)=w+b1*h0+r1*x0;

    for t=2:T
        h(t)=w+b1*h(t-1)+r1*x(t-1);
    end

    z=ret./sqrt(h);
    u=x-ks-ph*h-tau1*z-tau2*(z.^2-1);
    
end