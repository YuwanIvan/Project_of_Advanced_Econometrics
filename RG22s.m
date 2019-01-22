function [thetah, se, lrx, lr, pai, rho, ...
    rho_plus, rho_mins, h, z, u, Hmat, Jmat] = RG22s(ret,x, print)
    % Estimating the full model, that is, RealGARCH(2,2) with lagged returns added.
    % Output:
    % [thetah, se, lrx, lr, pai, rho, rho_plus, rho_mins, h, z, u, Hmat, Jmat]
    
    T=length(ret);
    options=optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0=[0; 1.45; -0.46; 0.42; -0.41; -0.18; 1.03; -0.07; 0.07; 0.15; 0];
    
    A=[0,0,0,0,0,0,0,0,0,-1,0];
    B=0;
    
    % minimise the negative log-likelihood function
    [thetah,~,~,~,~,~,hess]=fmincon(@(theta) -RG22s_logL(theta,ret,x),...
        theta0,A,B,[],[],[],[],[],options);
    Hmat=hess/T;
    Jmat=RG22s_J(thetah,ret,x);
    var=inv(Hmat)*Jmat*inv(Hmat)/T;
    se=diag(var);
    lrx=RG22s_logL(thetah,ret,x);
    lr=RG22s_logLp(thetah,ret,x);
    pai=thetah(11)+thetah(2)+thetah(3)+thetah(7)*(thetah(4)+thetah(5));
    [h,z,u]=RG22s_out(thetah,ret,x);
    item=thetah(8)*z+thetah(9)*(z.^2-1)+u;
    v0=corrcoef(item,z);
    rho=v0(2,1);
    vplus=corrcoef(item((z>0),1),z((z>0),1));
    rho_plus=vplus(2,1);
    vmins=corrcoef(item((z<0),1),z((z<0),1));
    rho_mins=vmins(2,1);
    
    if print==1
        disp('RealGARCH(2,2) Model with lagged returns');
        disp(' ');
        disp('[omega, beta1, beta2, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u, alpha] =');
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
        disp('rho =');
        disp(rho);
        disp(' ');
        disp('rho+ =');
        disp(rho_plus);
        disp(' ');
        disp('rho- =');
        disp(rho_mins);  
        
        disp(' ');
    end
    

    
end