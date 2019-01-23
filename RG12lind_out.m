function [h,z,u] = RG12lind_out(theta,ret,x)
    

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


end