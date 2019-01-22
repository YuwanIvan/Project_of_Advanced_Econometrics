function logL=G11_logL(theta,ret)
    % This is the log-likelihood function of 
    % GARCH(1,1) with lagged returns added.
    % The output is a scalar.
    
    logL=sum(G11_logf(theta,ret));
    
end