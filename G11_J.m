function J=G11_J(theta0,ret)
    % This is the J (outer-product) matrix of 
    % GARCH(1,1) with lagged returns added.
    % The output is a K by K matrix.
    
    T=length(ret);
    K=length(theta0);
    jt=numgrad(@G11_logf,theta0,ret);
    Jt=zeros(K,K);
    for t=1:T
        Jt=Jt+jt(t,:)'*jt(t,:);
    end      
    J=Jt/T;
end