function [z_sim, vz] = RG12N(ret,x)
    % News Impact Curve from RealGARCH(1,2)
    % Output:
    % [z_sim, vz]
    
    T = length(ret);
    options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');
    
    % starting values
    theta0 = [0.06; 0.7; 0.42; -0.18; -0.25; 1.03; 0.05; 0.085; 0.005; 0; 0.15];
    
    A = [0,0,0,0,0,0,0,0,0,0,-1];
    B = 0;
    
    % minimise the negative log-likelihood function
    thetah = fmincon(@(theta) -RG12N_logL(theta,ret,x),theta0,A,B,[],[],[],[],[],options);
    
    r1 = thetah(3);
    tau1 = thetah(7);
    tau2 = thetah(8);
    tau3 = thetah(9);
    tau4 = thetah(10);
    
    z_sim = -2:0.05:2;
    z_sim = z_sim';
    vz = (tau1 * z_sim + tau2 * (z_sim.^2-1) + tau3 * (z_sim.^3-3 * z_sim) + tau4 * (z_sim.^4-6 * z_sim.^2 + 3)) * r1;
    
    
    
end