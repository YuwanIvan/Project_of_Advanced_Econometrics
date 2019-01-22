function logL=RG11_logL(theta,ret,x)
    % RealGARCH(1,1) model.
    % The output is a scalar.
    
    logL=sum(RG11_logf(theta,ret,x));
    
end