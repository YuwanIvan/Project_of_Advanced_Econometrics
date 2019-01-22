function logL=RG12d_logL(theta,ret,x)
    % RealGARCH(1,2) model without leverage.
    % The output is a scalar.
    
    logL=sum(RG12d_logf(theta,ret,x));
    
end