function logL=RG22s_logL(theta,ret,x)
    % This is the log-likelihood function of the full
    % model, that is, RealGARCH(2,2) with lagged returns added.
    % The output is a scalar.
    
    logL=sum(RG22s_logf(theta,ret,x));
    
end