function logf = RG12lind_logf(theta,ret,x)
    % This is the observation-specific log-likelihood function of
    % RealGARCH(1,2) model without leverage.
    % The output is a T-dimensional vector.
    
    % getdata
    T = length(ret);
    w = theta(1);
    b1 = theta(2);
    r1 = theta(3);
    r2 = theta(4);
    ks = theta(5);
    ph = theta(6);
    sig2 = theta(7);

    
    % initial values
    h0 = var(ret);
    x00 = mean(x);
    x0 = mean(x);
    h = zeros(T,1);
    h(1) = w + b1 * (h0) + r1 * (x0) + r2 * (x00);
    h(2) = w + b1 * (h(1)) + r1 * (x(1)) + r2 * (x0);


    % iteration
    for t = 3:T
        h(t) = w + b1 * (h(t-1)) + r1 * (x(t-1)) + r2 * (x(t-2));
    end


    % output
    z = ret ./ sqrt(h);
    u = x - ks - ph * (h);
    logf =  - 0.5 * log(2 * pi) - 0.5 * log(h) - (ret.^2) ./ (2 * h) - 0.5 * log(2 * pi * sig2) - 0.5*(u.^2) ./ sig2;
    
end