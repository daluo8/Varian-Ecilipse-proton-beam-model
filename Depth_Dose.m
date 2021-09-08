format long
close all;
%% Initialize
E0=200;%MeV
% alpha=2.623*10^-3;
% p=1.735;
% Rcsda=alpha*E0^p;
% b1=15.14450027;b2=29.8440076;
% g1=0.001260021;g2=0.003260031;
% Rcsda=6.8469e-3*E0*(1+b1-b1*exp(-g1*E0)+b2-b2*exp(-g2*E0));
a1=6.94656e-3;a2=8.13116e-4;a3=-1.21068e-6;a4=1.053e-9;
Rcsda=a1*E0+a2*E0^2+a3*E0^3+a4*E0^4;
rho=1;
tau_ini=0.;
x=0:0.1:45;
Ipp=Depth_Dosepp(E0,Rcsda,tau_ini,x);
Isec=Depth_Dosesec(E0,Rcsda,tau_ini,x);
[Irp,Iheavy]=Depth_Doserp(E0,Rcsda,tau_ini,x);
I=Ipp+Isec+Irp+Iheavy;

Nabs=rho/E0*sum(I);
Ipp=Ipp/Nabs;
I=Ipp+Isec+Irp+Iheavy;
% I=I./max(I);
% I=I./max(I);
figure
plot(x,I,'-')
hold on
plot(x,Ipp+Irp,'-')
plot(x,Isec,'-')
plot(x,Iheavy,'-')
legend('total','primary','secondary','recoil','Location','NorthWest')
xlabel('depth[cm]')
ylabel('dose')
%% plot MC sim curve
m=csvread('200MeV0%.csv',8,2);
x1=m(:,1)./100;
y1=m(:,2);
y1=y1./max(y1);
y1=flipud(y1);
figure
plot(x1,y1,'--','Color',[0 1 0])
hold on
plot(x,I./max(I),'-')
