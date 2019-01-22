function logL=RG22d_logL(theta,ret,x)
 
    % RealGARCH(2,2) without leverage..
    % The output is a scalar.
    
    logL=sum(RG22d_logf(theta,ret,x));
    
end