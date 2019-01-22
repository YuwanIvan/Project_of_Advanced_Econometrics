function [skew, kurt]=CumRet(theta,N)
    % Calculate the skewness and kurtosis of cumulative returns by
    % simulation
    % The sample size of simulation is T
    % the out put is a vector of skewness and kurtosis
    
    T=500;
    w=theta(1);
    b1=theta(2);

    r1=theta(3);
    r2=theta(4);
    ks=theta(5);
    ph=theta(6);
    tau1=theta(7);
    tau2=theta(8);
    sig2=theta(9);

    skw=zeros(251,N);
    kur=zeros(251,N);
    
    for n=1:1:N
        % initial values

        h0=0.883;
        x00=0.798;
        x0=0.798;


        h=zeros(T,1);
        x=zeros(T,1);

        z=randn(T,1);
        u=sqrt(sig2)*randn(T,1);

        h(1)=exp(w+b1*log(h0)+r1*log(x0)+r2*log(x00));
        x(1)=exp(ks+ph*log(h(1))+tau1*z(1)+tau2*(z(1)^2-1)+u(1));
        h(2)=exp(w+b1*log(h(1))+r1*log(x(1))+r2*log(x0));
        x(2)=exp(ks+ph*log(h(2))+tau1*z(2)+tau2*(z(2)^2-1)+u(2));

        for t=3:T
            h(t)=exp(w+b1*log(h(t-1))+r1*log(x(t-1))+...
                r2*log(x(t-2)));
            x(t)=exp(ks+ph*log(h(t))+tau1*z(t)+tau2*(z(t)^2-1)+u(t));
        end

        ret=z.*sqrt(h);

        skw(1,n)=skewness(ret);
        kur(1,n)=kurtosis(ret);
        
        for K=1:1:250
            temp=zeros(T-K,1);
            for t=1:1:(T-K)
                temp(t)=sum(ret(t:t+K,1));
            end
            skw(K+1,n)=skewness(temp);
            kur(K+1,n)=kurtosis(temp);
        end
    end
    
    skew=zeros(251,1);
    kurt=zeros(251,1);
    for i=1:251
        skew(i)=mean(skw(i,:));
        kurt(i)=mean(kur(i,:));
    end
    
    
end