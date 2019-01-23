function logL=RG12N_logL(theta,ret,x)
    % RealGARCH(1,2) model with 4th leverage.
    % The output is a scalar.
    
    logL = sum(RG12N_logf(theta,ret,x));
    
end