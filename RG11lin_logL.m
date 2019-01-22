function logL=RG11lin_logL(theta,ret,x)
    % RealGARCH(1,1) model with linear specification.
    % The output is a scalar.
    
    logL=sum(RG11lin_logf(theta,ret,x));
    
end