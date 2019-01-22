function [insamp, outsamp, outpart]=LRtest(retin,xin,retout,xout)
    % log-likelihood ratio test
    % sytax:
    % [insamp, outsamp, outpart]=LRtest(retin,xin,retout,xout)
    
    tg=G11(retin,0);
    t11=RG11(retin,xin,0);
    t12=RG12(retin,xin,0);
    t21=RG21(retin,xin,0);
    t22=RG22(retin,xin,0);
    t22d=RG22d(retin,xin,0);
    t22s=RG22s(retin,xin,0);
    
    thetag=[tg(1);tg(2);zeros(8,1);tg(3)];
    theta11=[t11(1);t11(2);0;t11(3);0;t11(4);t11(5:8,1);0];
    theta12=[t12(1);t12(2);0;t12(3:9,1);0];
    theta21=[t21(1:4,1);0;t21(5:9,1);0];
    theta22=[t22;0];
    theta22d=[t22d(1:7,1);0;0;t22d(8);0];
    theta22s=t22s;
    
    in1=RG22s_logL(theta11,retin,xin);
    in2=RG22s_logL(theta12,retin,xin);
    in3=RG22s_logL(theta21,retin,xin);
    in4=RG22s_logL(theta22,retin,xin);
    in5=RG22s_logL(theta22d,retin,xin);
    in6=RG22s_logL(theta22s,retin,xin);
    
    out1=RG22s_logL(theta11,retout,xout);
    out2=RG22s_logL(theta12,retout,xout);
    out3=RG22s_logL(theta21,retout,xout);
    out4=RG22s_logL(theta22,retout,xout);
    out5=RG22s_logL(theta22d,retout,xout);
    out6=RG22s_logL(theta22s,retout,xout);
    
    pg=RG22s_logLp(thetag,retout,xout);
    p1=RG22s_logLp(theta11,retout,xout);
    p2=RG22s_logLp(theta12,retout,xout);
    p3=RG22s_logLp(theta21,retout,xout);
    p4=RG22s_logLp(theta22,retout,xout);
    p5=RG22s_logLp(theta22d,retout,xout);
    p6=RG22s_logLp(theta22s,retout,xout);
    
    insamp=zeros(6,1);
    insamp(1)=in6-in1;
    insamp(2)=in6-in2;
    insamp(3)=in6-in3;
    insamp(4)=in6-in4;
    insamp(5)=in6-in5;
    insamp(6)=in6-in6;
    insamp=2*insamp;
    
    outsamp=zeros(6,1);
    outsamp(1)=out4-out1;
    outsamp(2)=out4-out2;
    outsamp(3)=out4-out3;
    outsamp(4)=out4-out4;
    outsamp(5)=out4-out5;
    outsamp(6)=out4-out6;
    outsamp=sqrt(length(retin)/length(retout))*outsamp;
    
    partmax=max([pg,p1,p2,p3,p4,p5,p6]');
    outpart=zeros(7,1);
    outpart(1)=partmax-pg;
    outpart(2)=partmax-p1;
    outpart(3)=partmax-p2;
    outpart(4)=partmax-p3;
    outpart(5)=partmax-p4;
    outpart(6)=partmax-p5;
    outpart(7)=partmax-p6;
    outpart=2*outpart;
    
    disp(' ');
    disp('In-sample: 11, 12, 21, 22, 22+, 22*');
    disp(insamp');
    disp(' ');
    disp('Out-of-sample: 11, 12, 21, 22, 22+, 22*');
    disp(outsamp');
    disp(' ');
    disp('Out-of-sample partial:G11, 11, 12, 21, 22, 22+, 22*');    
    disp(outpart');
    disp(' ');

    
end