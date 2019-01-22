function [thetah, se, lrx, lr, pai, rho, ...
    rho_plus, rho_mins, h, z, u, Hmat, Jmat] = RG12(ret,x, print)
    % Estimating RealGARCH(1,2)
    % Output:
    % [thetah, se, lrx, lr, pai, rho, rho_plus, rho_mins, h, z, u, Hmat, Jmat]
    
    T=length(ret);
    options=optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0=[0.06; 0.7; 0.42; -0.18; -0.25; 1.03; -0.07; 0.07; 0.15];
    
    A=[0,0,0,0,0,0,0,0,-1];
    B=0;
    
    % minimise the negative log-likelihood function
    [thetah,~,~,~,~,~,hess]=fmincon(@(theta) -RG12_logL(theta,ret,x),...
        theta0,A,B,[],[],[],[],[],options);
    Hmat=hess/T;
    Jmat=RG12_J(thetah,ret,x);
    var=inv(Hmat)*Jmat*inv(Hmat)/T;
    se=diag(var);
    lrx=RG12_logL(thetah,ret,x);
    lr=RG12_logLp(thetah,ret,x);
    pai=thetah(2)+thetah(6)*(thetah(3)+thetah(4));
    [h,z,u]=RG12_out(thetah,ret,x);
    item=thetah(7)*z+thetah(8)*(z.^2-1)+u;
    v0=corrcoef(item,z);
    rho=v0(2,1);
    vplus=corrcoef(item((z>0),1),z((z>0),1));
    rho_plus=vplus(2,1);
    vmins=corrcoef(item((z<0),1),z((z<0),1));
    rho_mins=vmins(2,1);
    
    if print==1
        disp('RealGARCH(1,2) Model');
        disp(' ');
        disp('[omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u] =');
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