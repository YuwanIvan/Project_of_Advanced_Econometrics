function [thetah, se, lrx, lr, pai, h, z, u, Hmat, Jmat] = RG22d(ret,x, print)
    % Estimating the full model, that is, RealGARCH(2,2) without leverage.
    % Output:
    % [thetah, se, lrx, lr, pai, h, z, u, Hmat, Jmat]
    
    T=length(ret);
    options=optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0=[0; 1.45; -0.46; 0.42; -0.41; -0.18; 1.03;  0.15];
    
    A=[0,0,0,0,0,0,0,-1];
    B=0;
    
    % minimise the negative log-likelihood function
    [thetah,~,~,~,~,~,hess]=fmincon(@(theta) -RG22d_logL(theta,ret,x),...
        theta0,A,B,[],[],[],[],[],options);
    Hmat=hess/T;
    Jmat=RG22d_J(thetah,ret,x);
    var=inv(Hmat)*Jmat*inv(Hmat)/T;
    se=diag(var);
    lrx=RG22d_logL(thetah,ret,x);
    lr=RG22d_logLp(thetah,ret,x);
    pai=thetah(2)+thetah(3)+thetah(7)*(thetah(4)+thetah(5));
    [h,z,u]=RG22d_out(thetah,ret,x);
    
    if print==1
        disp('RealGARCH(2,2) Model without leverage');
        disp(' ');
        disp('[omega, beta1, beta2, gamma1, gamma2, xi, phi, sig2_u] =');
        disp(thetah');
        disp(' ');
        disp('std. err. =');
        disp(se');
        disp(' ');
        disp('l(r,x) =');
        disp(lrx);
        disp(' ');
        disp('l(r) =');
        disp(lr);
        disp(' ');
        disp('pi =');
        disp(pai);

        disp(' ');
    end
    

    
end