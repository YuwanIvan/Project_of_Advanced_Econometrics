function [thetah, se, lrx, lr, pai, rho, rho_plus, rho_mins, h, z, u, Hmat, Jmat] = RG12lin(ret, x, print)
    % Estimating RealGARCH(1,2) with linear specification
    % Output:
    % [thetah, se, lrx, lr, pai, rho, rho_plus, rho_mins, h, z, u, Hmat, Jmat]
    
    %% Estimate
    T = length(ret);
    options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    % starting values
    theta0 = [0.06; 0.7; 0.42; -0.18; -0.25; 1.03; -0.07; 0.07; 0.15];
    % constraints
    A = [0, 0, 0, 0, 0, 0, 0, 0, -1];
    B = 0;
    % minimise the negative log-likelihood function
    [thetah, ~, ~, ~, ~, ~, hess] = fmincon(@(theta) -RG12lin_logL(theta,ret,x), theta0, A, B, [], [], [], [], [], options);
    
    %% output
    Hmat = hess/T;
    Jmat = RG12lin_J(thetah,ret,x);
    var = inv(Hmat) * Jmat * inv(Hmat) / T;
    se = diag(var);
    lrx = RG12lin_logL(thetah,ret,x);
    lr = RG12lin_logLp(thetah,ret,x);
    pai = thetah(2) + thetah(6) * (thetah(3)+thetah(4));
    [h,z,u] = RG12lin_out(thetah,ret,x);
    item = thetah(7) * z +
    thetah(8) * (z.^2 - 1) + u;
    v0 = corrcoef(item,z);
    rho = v0(2,1);
    vplus = corrcoef(item((z>0),1),z((z>0),1));
    rho_plus = vplus(2,1);
    vmins = corrcoef(item((z<0),1),z((z<0),1));
    rho_mins = vmins(2,1);
    
    if print==1
        disp('RealGARCH(1,2) Model with linear specification');
        disp(' ');
        disp('[omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u] =');
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