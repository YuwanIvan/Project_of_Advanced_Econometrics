function logL=RG12_logL(theta,ret,x)
    % RealGARCH(1,2) model.
    % The output is a scalar.
    
    logL=sum(RG12_logf(theta,ret,x));
    
end