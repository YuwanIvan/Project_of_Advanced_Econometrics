%*************************************************************************
%
% Assignment 3 for Advanced Econometrics (1)
%
%*************************************************************************


clear all;
close all;
clc;

delete assignment03.txt
diary assignment03.txt

fig=figure;
fig.PaperSize=[6,4.5];
fig.PaperPosition=[0.1, 0.1, 5.8, 4.3];

RandStream.setGlobalStream(RandStream('mt19937ar','seed',78364))

disp('***********************************************************************');
disp(' ');
disp('Assignment 3 for Advanced Econometrics (1)');
disp('Research Replication');
disp(' ');
disp('Original Work:');
disp('Hansen, P. R., Huang, Z., & Shek, H. H. (2012). Realized GARCH:');
disp('A joint model for returns and realized measures of volatility.');
disp('Journal of Applied Econometrics, 27(6), 877-906.');
disp(' ');
disp('***********************************************************************');
disp(' ');

%**************************************************************************
% Table 2
%**************************************************************************

clear all
disp('=================================================================');
disp(' ');
disp('Table 2: Log-linear Models for SPY');
disp(' ');

data=csvread('hhs_data\SPY.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);
ccret=data((date<20080000),4);

disp('(1) Open-to-close returns ');
disp(' ');

disp('---------------------');
G11(ocret,1);
disp('---------------------');
RG11(ocret,x,1);
disp('---------------------');
RG12(ocret,x,1);
disp('---------------------');
RG21(ocret,x,1);
disp('---------------------');
RG22(ocret,x,1);
disp('---------------------');
RG22d(ocret,x,1);
disp('---------------------');
RG22s(ocret,x,1);
disp('---------------------');

disp('(2) Close-to-close returns ');
disp(' ');

disp('---------------------');
G11(ccret,1);
disp('---------------------');
RG11(ccret,x,1);
disp('---------------------');
RG12(ccret,x,1);
disp('---------------------');
RG21(ccret,x,1);
disp('---------------------');
RG22(ccret,x,1);
disp('---------------------');
RG22d(ccret,x,1);
disp('---------------------');
RG22s(ccret,x,1);
disp('---------------------');


%**************************************************************************
% Table 3
%**************************************************************************
clear all;

disp('=================================================================');
disp(' ');
disp('Table 3: Log-linear RealGARCH(1,2) Model');
disp(' ');

Data_names = {'AA', 'AIG', 'AXP', 'BA', 'BAC', 'C', 'CAT', 'CVX', 'DD', 'DIS', 'GE', 'GM', 'HD', 'IBM', 'INTC', 'JNJ', 'JPM', 'KO', 'MCD', 'MMM', 'MRK', 'MSFT', 'PG',  'T', 'UTX', 'VZ', 'WMT', 'XOM', 'SPY'};
data_num = length(Data_names);

for data_id = 1 : data_num
	clear data date ocret x;

	data_name = Data_names{data_id};
	data_path = ['hhs_data','\' data_name, '.csv' ];

	disp('----------------------------');
	disp(data_name);
	data=csvread(data_path);
	date=data(:,1);
	ocret=data((date<20080000),2);
	x=data((date<20080000),3);
	[~,~,~,~,~,~,~,~,~,~,~,~,~]=RG12(ocret,x,1);	
end

%**************************************************************************
% Figure 1
%**************************************************************************
clear all

fig=figure;
fig.PaperSize=[6,4.5];
fig.PaperPosition=[0.1, 0.1, 5.8, 4.3];

data=csvread('hhs_data/SPY.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);
[z,vz1]=RG12N(ocret,x);

data=csvread('hhs_data/IBM.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);
[z,vz2]=RG12N(ocret,x);

plot(z,vz1,'-',z,vz2,'--');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
title('Figure 1: News Impact Curve');
xlabel('z');
ylabel('v(z)');
legend('SPY','IBM','Location','northeast');
print('-dpdf','Figure_1.pdf');


%**************************************************************************
% Table 4
%**************************************************************************

clear all;

disp('=================================================================');
disp(' ');
disp('Table 4: In-sample and out-of-sample likelihood ratio test');
disp(' ');

Data_names = {'AA', 'AIG', 'AXP', 'BA', 'BAC', 'C', 'CAT', 'CVX', 'DD', 'DIS', 'GE', 'GM', 'HD', 'IBM', 'INTC', 'JNJ', 'JPM', 'KO', 'MCD', 'MMM', 'MRK', 'MSFT', 'PG',  'T', 'UTX', 'VZ', 'WMT', 'XOM', 'SPY'};
data_num = length(Data_names);


for data_id = 1 : data_num
	clear data date ocret ocret_out x x_out;

	data_name = Data_names{data_id};
    data_path = ['hhs_data','\' data_name, '.csv' ];

	disp('----------------------------');
	disp(data_name);
	data=csvread(data_path);
	date=data(:,1);
	ocret=data((date<20080000),2);
	ocret_out=data((date>20080000),2);
	x=data((date<20080000),3);
	x_out=data((date>20080000),3);
	[~,~,~]=LRtest(ocret,x,ocret_out,x_out);	
end





%**************************************************************************
% Figure 2
%**************************************************************************

clear all

data=csvread('hhs_data/SPY.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);

[~,~,~,~,~,~,~,~,h1]=RG12lin(ocret,x,0);
[~,~,~,~,~,~,~,~,h2]=RG12(ocret,x,0);

figure
title('Figure 2: Heteroskedasticity in measurement equation');

subplot(1,2,1);
plot(h1,x,'r.');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
set(gcf,'position',[200,100,800,350]);
title('Figure 2(1): Linear Model: x against h');
refline([0,2]);
refline([0,4]);
refline([0,6]);
refline([0,8]);
refline([0,10]);
refline([0,12]);
refline([0,14]);
xlabel('h');
ylabel('x');
print('-dpdf','Figure_2(1).pdf');

subplot(1,2,2);
plot(log(h2),log(x),'r.');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
set(gcf,'position',[200,100,800,350]);
title('Figure 2(2): Log-linear Model: log(x) against log(h)');
refline([0,1]);
refline([0,2]);
refline([0,3]);
refline([0,-1]);
refline([0,-2]);
refline([0,-3]);
xlabel('log(h)');
ylabel('log(x)');
print('-dpdf','Figure_2(2).pdf');



%**************************************************************************
% Figure 3
%**************************************************************************

clear all

fig=figure;
fig.PaperSize=[6,4.5];
fig.PaperPosition=[0.1, 0.1, 5.8, 4.3];

data=csvread('hhs_data/SPY.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);

[~,~,~,~,~,~,~,~,~,z4,u4]=RG12(ocret,x,0);
[~,~,~,~,~,~,~,~,~,z2,u2]=RG12lin(ocret,x,0);
[~,z1,u1]=RG12lind(ocret,x);
[~,z3,u3]=RG12d(ocret,x);

subplot(2,2,1);
plot(z1,u1,'b.');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
set(gca,'xtick',[],'xticklabel',[]);
set(gca,'ytick',[],'yticklabel',[]);
set(gcf,'position',[200,100,800,800]);
title('Figure 3(1):Linear specification without leverage function');
xlabel('z');
ylabel('u');
print('-dpdf','Figure_3(1).pdf');

subplot(2,2,2);
plot(z2,u2,'b.');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
set(gca,'xtick',[],'xticklabel',[]);
set(gca,'ytick',[],'yticklabel',[]);
set(gcf,'position',[200,100,800,800]);
title('Figure 3(2):Linear specification with quadratic leverage function');
xlabel('z');
ylabel('u');
print('-dpdf','Figure_3(2).pdf');

subplot(2,2,3);
plot(z3,u3,'b.');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
set(gca,'xtick',[],'xticklabel',[]);
set(gca,'ytick',[],'yticklabel',[]);
set(gcf,'position',[200,100,800,800]);
title('Figure 3(3):Log-linear specification without leverage function');
xlabel('z');
ylabel('u');
print('-dpdf','Figure_3(3).pdf');

subplot(2,2,4);
plot(z4,u4,'b.');
set(gca,'XAxisLocation','origin');
set(gca,'YAxisLocation','origin');
set(gca,'xtick',[],'xticklabel',[]);
set(gca,'ytick',[],'yticklabel',[]);
set(gcf,'position',[200,100,800,800]);
title('Figure 3(4):Log-linear specification with quadratic leverage function');
xlabel('z');
ylabel('u');
print('-dpdf','Figure_3(4).pdf');

%**************************************************************************
% Table 5
%**************************************************************************

clear all;

data=csvread('hhs_data\SPY.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);



[~,se,~,~,~,~,~,~,~,~,~,H,J]=RG12lin(ocret,x,0);
T=length(ocret);
se_I=sqrt(diag(inv(H)/T));
se_J=sqrt(diag(inv(J)/T));

disp('=================================================================');
disp(' ');
disp('Table 5: Standard errors for RealGARCH(1,2)');
disp(' ');
disp('(1) Linear Model');
disp(' ');
disp('std. err. based on I');
disp('omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u');
disp(se_I');
disp(' ');
disp('std. err. based on J');
disp('omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u');
disp(se_J');
disp(' ');
disp('std. err. based on sandwich');
disp('omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u');
disp(se');
disp(' ');
disp(' ');


[~,se,~,~,~,~,~,~,~,~,~,H,J]=RG12(ocret,x,0);
T=length(ocret);
se_I=sqrt(diag(inv(H)/T));
se_J=sqrt(diag(inv(J)/T));

disp('(2) Log-linear Model');
disp(' ');
disp('std. err. based on I');
disp('omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u');
disp(se_I');
disp(' ');
disp('std. err. based on J');
disp('omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u');
disp(se_J');
disp(' ');
disp('std. err. based on sandwich');
disp('omega, beta1, gamma1, gamma2, xi, phi, tau1, tau2, sig2_u');
disp(se');
disp(' ');
disp(' ');


%**************************************************************************
% Figure 4
%**************************************************************************

clear all;

fig=figure;
fig.PaperSize=[6,4.5];
fig.PaperPosition=[0.1, 0.1, 5.8, 4.3];



RandStream.setGlobalStream(RandStream('mt19937ar','seed',78364))

data=csvread('hhs_data/SPY.csv');
date=data(:,1);
ocret=data((date<20080000),2);
x=data((date<20080000),3);

thetah=RG12(ocret,x,0);
[skew, kurt]=CumRet(thetah,1000);
lag=0:1:250;

subplot(2,1,1);
plot(lag',skew,'b-');
set(gca,'XAxisLocation','origin');
set (gca, 'XTick', 0:25:250);
set(gca,'YAxisLocation','origin');
set (gca, 'YTick', -0.30:0.05:0);
title('Figure 4(1): Skewness');
set(gcf,'position',[200,100,400,900]);
refline([0,-0.05]);
refline([0,-0.10]);
refline([0,-0.15]);
refline([0,-0.20]);
refline([0,-0.25]);
refline([0,-0.30]);
xlabel('lags');
ylabel('skewness');
print('-dpdf','Figure_4(1).pdf');

subplot(2,1,2);
plot(lag',kurt,'b-');
ylim([2,5.5]);
set (gca, 'XTick', 0:25:250);
set(gca,'YAxisLocation','origin');
set (gca, 'YTick', 2:0.5:5.50);
set(gcf,'position',[200,100,400,900]);
refline([0,2.00]);
refline([0,2.50]);
refline([0,3.00]);
refline([0,3.50]);
refline([0,4.00]);
refline([0,4.50]);
refline([0,5.00]);
refline([0,5.50]);
title('Figure 4(2): Kurtosis');
xlabel('lags');
ylabel('kurtosis');
print('-dpdf','Figure_4(2).pdf');

diary off

clear all;
close all;

% The end
