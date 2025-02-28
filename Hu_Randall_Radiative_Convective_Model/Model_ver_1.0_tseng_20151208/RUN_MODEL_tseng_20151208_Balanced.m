
clear all;

N=10;
yy  =linspace(1,10,N);
xx  =linspace(298,304,N);
%sm  =linspace(302*1004,300*1004,N);
%sbt =linspace(302*1004,298*1004,N);
%qm  =linspace(0.008,0.024,N);
%qbt =linspace(0.004,0.020,N);
%mc  =linspace(0.001,0,N);



[SST,Vs]=meshgrid(xx,yy);


c=0.6;
SM  = zeros(size(Vs,2),size(SST,2));
SBt = zeros(size(Vs,2),size(SST,2));
qM  = zeros(size(Vs,2),size(SST,2));
qBt = zeros(size(Vs,2),size(SST,2));
PR  = zeros(size(Vs,2),size(SST,2));




X0(1)=3.1*10^5;  %first guess of SM
X0(2)=3*10^5;  %first guess of SBt
X0(3)=0.02;    %first guess of qM
X0(4)=0.02;    %first guess of qBt
X0(5)=1e-4;       %first guess of Mc


%%

for j=1:size(Vs,2)
    
    %X0(1)=sm(j);      %first guess of SM
    %X0(2)=sbt(j);      %first guess of SBt
    %X0(3)=qm(j);      %first guess of qM
    %X0(4)=qbt(j);      %first guess of qBt
    %X0(5)=mc(j);      %first guess of 
    
    for i=1:size(SST,2)

    fun2=@(x) balanced(x,Vs(i,j),SST(i,j));
    x=fsolve(fun2,X0);


    SM (i,j) = x(1);
    SBt(i,j) = x(2); 
    qM (i,j) = x(3); qM(qM<0)=0;
    qBt(i,j) = x(4); qBt(qBt<0)=0;
    PR (i,j) = c*x(3)*x(5)*3600;
    end
end

%%
figure(1)
contourf(SST,Vs,SBt/1004);
colorbar
title('Mixed Layer Temperature')
xlabel('SST')
ylabel('surface wind')
set(gca,'fontsize',15)

figure(2)
contourf(SST,Vs,PR);
colorbar
title('Precipitation Rate (mm/hr)')
xlabel('SST')
ylabel('surface wind')
set(gca,'fontsize',15)

