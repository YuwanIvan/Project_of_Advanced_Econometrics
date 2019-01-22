function J=RG22d_J(theta0,ret,x)
    % RealGARCH(2,2) without leverage.
    % The output is a K by K matrix.
    
    T=length(ret);
    K=length(theta0);
    jt=numgrad(@RG22d_logf,theta0,ret,x);
    Jt=zeros(K,K);
    for t=1:T
        Jt=Jt+jt(t,:)'*jt(t,:);
    end      
    J=Jt/T;
end