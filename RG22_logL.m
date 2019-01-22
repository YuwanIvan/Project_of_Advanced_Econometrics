function logL=RG22_logL(theta,ret,x)
    % RealGARCH(2,2) model.
    % The output is a scalar.
    
    logL=sum(RG22_logf(theta,ret,x));
    
end