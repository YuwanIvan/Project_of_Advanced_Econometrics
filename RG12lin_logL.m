function logL=RG12lin_logL(theta,ret,x)
    % RealGARCH(1,2) model with linear specification.
    % The output is a scalar.
    
    logL=sum(RG12lin_logf(theta,ret,x));
    
end