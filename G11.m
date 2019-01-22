function [thetah, se, lr, pai, h, z, Hmat, Jmat] = G11(ret,print)
    % Estimating GARCH(1,1) with lagged returns added.
    % Output:
    % [thetah, se, lr, pai, h, z, Hmat, Jmat]
    
    T=length(ret);
    options=optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0=[0.05; 0.96; 0.03];
    
    A=[0,1,1; 0,-1,0; 0,0,-1];
    B=[1;0;0];
   
    % minimise the negative log-likelihood function
    [thetah,~,~,~,~,~,hess]=fmincon(@(theta) -G11_logL(theta,ret),...
        theta0,A,B,[],[],[],[],[],options);
    Hmat=hess/T;
    Jmat=G11_J(thetah,ret);
    var=inv(Hmat)*Jmat*inv(Hmat)/T;
    se=diag(var);
    lr=G11_logL(thetah,ret);
    pai=thetah(2)+thetah(3);
    [h,z]=G11_out(thetah,ret);

    
    if print==1
        disp('GARCH(1,1) Model');
        disp(' ');
        disp('[omega, beta1, alpha] =');
        disp(thetah');
        disp(' ');
        disp('std. err. =');
        disp(se');
        disp(' ');
        disp('l(r) =');
        disp(lr);
        disp(' ');
        disp('pi =');
        disp(pai);
        disp(' ');
    end
    

    
end