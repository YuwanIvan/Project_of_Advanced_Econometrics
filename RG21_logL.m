function logL=RG21_logL(theta,ret,x)
    % RealGARCH(2,1) model.
    % The output is a scalar.
    
    logL=sum(RG21_logf(theta,ret,x));
    
end