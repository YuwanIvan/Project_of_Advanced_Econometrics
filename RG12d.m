function [h, z, u] = RG12d(ret,x)
    % Estimating RealGARCH(1,2) without leverage
    % Output:
    % [h, z, u]
    
    T=length(ret);
    options=optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0=[0.06; 0.7; 0.42; -0.18; -0.25; 1.03;  0.15];
    
    A=[0,0,0,0,0,0,-1];
    B=0;
    
    % minimise the negative log-likelihood function
    thetah=fmincon(@(theta) -RG12d_logL(theta,ret,x),...
        theta0,A,B,[],[],[],[],[],options);
    [h,z,u]=RG12d_out(thetah,ret,x);
        
   

    
end