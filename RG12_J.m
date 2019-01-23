function J=RG12_J(theta0,ret,x)
    % This is the J (outer-product) matrix of 
    % RealGARCH(1,2)
    % The output is a K by K matrix.
    
    T = length(ret);
    K = length(theta0);
    jt = numgrad(@RG12_logf,theta0,ret,x);
    Jt = zeros(K,K);
    
    for t = 1:T
        Jt = Jt + jt(t,:)' * jt(t,:);
    end

    J=Jt/T;
end