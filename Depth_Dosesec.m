function Isec = Depth_Dosesec(E0,Rcsda,tau_ini,x)
format long
% tau_ini=0;
tau_straggle=0.0127*Rcsda^0.9352;
switch E0
    case E0<7
        tau_heavy=0.2e-6;
    case E0>7&&E0<=20.12
        tau_heavy=0.2e-6+1.10822*(E0-7)/(20.12-7);
    otherwise
        tau_heavy=1.108222-0.001170874*(E0-20.12);
end
tau=sqrt(tau_ini^2+tau_straggle^2+tau_heavy^2);
%% get phi array
% phi=Fraction_of_primary(E0,Rcsda);
%% calculate Ci
C1=5.7088-0.0046297*E0;
C2=0.4862+0.0014*E0;
C3=2.0796-0.0020806*E0;
C4=0.8106-0.0014*E0;
k1=9.8492+0.0017483*E0;
%select the value of epi[mm] according to E0
switch E0
    case E0<50
        epi=2*10^-6;
    case E0>=50&&E0<=70
        epi=2*10^-6+9.9*10^-6*(E0-50);
    otherwise
        epi=2*10^-4+3.4574*10^-4*(E0-70)+1.22*10^(-4)*(E0-70)^2;
end
zmax=Rcsda+epi/10;
%% calculate Ii
% x=0:0.1:40;
if E0<=20.12%MeV
    xshift=0.255*exp(-2*pi^2*(20.12-E0)^2/20.12^2);
else
    xshift=0.255*exp(-(E0-20.12)^2/106.87541^2);
end
% xshift=0;
xs=x+xshift;
I1=(C1*tau_straggle-C4*(tau/Rcsda)^2*(Rcsda+xs)).*1./(sqrt(2*pi)*tau).*exp(-(Rcsda-xs).^2/(2*tau^2));
I2=(C2+C4/sqrt(pi)*(tau/Rcsda)^2)*1/2*(1+erf((Rcsda-xs)/(sqrt(2)*tau)));
I3=C3*exp(2*((k1*tau_ini)/zmax)^2-2*k1*(Rcsda-xs)/zmax)*1/2.*(1+erf((Rcsda-xs)/(sqrt(2)*tau)-sqrt(2)*k1*tau/zmax));
I4=C4*(xs/Rcsda).^2*1/2.*(1+erf((Rcsda-xs)/(sqrt(2)*tau)));

%% calculate Clani,Rlani and taulani
% in Landau terms we use z1=z-znet and Rcsda_patient=Rcsda-znet
xnet=0;%water-equivalent thickness of material in nozzle
Rcsda_patient=Rcsda-xnet;
Clan1=-6.850*10^-6*Rcsda_patient^3+4.779*10^-4*Rcsda_patient^2 ...
      -8.199*10^-3*Rcsda_patient-2.545*10^-2;
switch E0%here assump Epatient=E0
    case E0<120
        Clan2=4.3582*10^-5*E0;
    case E0>=120&&E0<=168
        Clan2=2.8476*10^-3*(E0-118)^0.877;
    otherwise
        Clan2=0.088;
end
if E0<68
    Rlan1=0.7*Rcsda_patient;
else
    Rlan1=(0.81209-0.0016484*E0)*Rcsda_patient;
end
Rlan2=(3.19+0.00161*E0)*(1-exp(-(E0/165.8)^2));
% taulan1=0.7071*Rlan1+0.0492*tau_ini;
taulan1=0.7071*Rlan1+0.0492*tau_ini+tau_heavy;
taulan2=2.4*tau;
k2=2.56;
%% calculate Lani
x1=xs-xnet;
Ilan1=Clan1*(erf(2*x1/Rlan1)+erf(k2*(Rlan1-x1)/(sqrt(2)*taulan1)));
Ilan2=Clan2/Rlan2^2*((Rlan2^2-x1.^2-taulan2^2/sqrt(pi))*1/2.*(1+erf((Rlan2-x1)/(sqrt(2)*taulan2))) ...
    +taulan2*(x1+Rlan2)/sqrt(2*pi).*exp(-(Rlan2-x1).^2/(2*taulan2^2)));
Ilan1=Ilan1;
Ilan2=-Ilan2;

% legend('Ilan1','Ilan2')
%% calculate I
Mpc2=938.272;
Eth=7;
xi=((E0-Eth)/Mpc2)^1.032;
Isec=0.958*xi*x/Rcsda.*(I1+I2+I3+I4+Ilan1+Ilan2);

end

