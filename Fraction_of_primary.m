function phi = Fraction_of_primary(E0,Rcsda)

Mpc2=938.272;%MeV
Eth=7;%MeV
xi=((E0-Eth)/Mpc2)^1.032;
tau=0.0127*Rcsda^0.9352;

x=0:0.1:40;
n=length(x);
phi=zeros(n,1);
phi=(1-xi*x./Rcsda)*1/2.*(1+erf((Rcsda-x)./(sqrt(2)*tau)));
end

