function [thetah, se, lrx, lr, pai, rho, ...
    rho_plus, rho_mins, h, z, u, Hmat, Jmat] = RG11lin(ret,x, print)
    % Estimating RealGARCH(1,1) with linear specification
    % Output:
    % [thetah, se, lrx, lr, pai, rho, rho_plus, rho_mins, h, z, u, Hmat, Jmat]
    
    T=length(ret);
    options=optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0=[0.18; 0.55; 0.42; -0.42; 1.03; -0.11; 0.04; 0.28];
    
    A=[0,0,0,0,0,0,0,-1];
    B=0;
    
    % minimise the negative log-likelihood function
    [thetah,~,~,~,~,~,hess]=fmincon(@(theta) -RG11lin_logL(theta,ret,x),...
        theta0,A,B,[],[],[],[],[],options);
    Hmat=hess/T;
    Jmat=RG11lin_J(thetah,ret,x);
    var=inv(Hmat)*Jmat*inv(Hmat)/T;
    se=diag(var);
    lrx=RG11lin_logL(thetah,ret,x);
    lr=RG11lin_logLp(thetah,ret,x);
    pai=thetah(2)+thetah(5)*(thetah(3));
    [h,z,u]=RG11lin_out(thetah,ret,x);
    item=thetah(6)*z+thetah(7)*(z.^2-1)+u;
    v0=corrcoef(item,z);
    rho=v0(2,1);
    vplus=corrcoef(item((z>0),1),z((z>0),1));
    rho_plus=vplus(2,1);
    vmins=corrcoef(item((z<0),1),z((z<0),1));
    rho_mins=vmins(2,1);
    
    if print==1
        disp('RealGARCH(1,1) Model with linear specification');
        disp(' ');
        disp('[omega, beta1, gamma1, xi, phi, tau1, tau2, sig2_u] =');
        disp(thetah');
        disp('std. err. =');
        disp(se');
        disp('l(r,x) =');
        disp(lrx);
        disp('l(r) =');
        disp(lr);
        disp('pi =');
        disp(pai);
        disp('rho =');
        disp(rho);
        disp('rho+ =');
        disp(rho_plus);
        disp('rho- =');
        disp(rho_mins);  
        disp(' ');
    end
    

    
end