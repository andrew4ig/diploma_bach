clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;

%initial paramets
a=3e-9;                 %size of chanel, 3nm
b=2e-9;                 %size of barrier, 2nm
U0=1*e;                 %height of barrier, 1 eV

eps1=12.90;
m1e=0.067*m0;           %eff. mass in GaAs - e, kg
m1h=0.082*m0;           %eff. mass in GaAs - holes, kg
Eg1=1.42*e;
Ev1=0;
Ec1=Ev1+Eg1;

eps2=12.90-2.84*0.5;
m2e=(0.063+0.083*0.5)*m0;%eff. mass in AlGaAs(0.5) - e, kg
m2h=(0.082+0.068*0.5)*m0;%eff. mass in AlGaAs(0.5) - holes, kg
Eg2=(1.9+0.125*0.5+0.143*0.5^2)*e;
Ev2=0;
Ec2=Ev2+Eg2;

T=300;                  %room temperature, K

L=3*a+2*b;
N=L*1e9*10;
x=linspace(-L/2,L/2,N);

mh=[m1h*ones(1,ceil(a*N/L)),m2h*ones(1,b/L*N),m1h*ones(1,ceil(a*N/L)),...
   m2h*ones(1,b/L*N),m1h*ones(1,ceil(a/L*N))];
me=[m1e*ones(1,ceil(a*N/L)),m2e*ones(1,b/L*N),m1e*ones(1,ceil(a*N/L)),...
   m2e*ones(1,b/L*N),m1e*ones(1,ceil(a*N/L))];

phi=0;
eps=[eps1*ones(1,ceil(a*N/L)),eps2*ones(1,b/L*N),eps1*ones(1,ceil(a*N/L)),...
   eps2*ones(1,b/L*N),eps1*ones(1,ceil(a/L*N))];
% Ec=[Ec1*ones(1,ceil(a*N/L)),Ec2*ones(1,b/L*N),Ec1*ones(1,ceil(a*N/L)),...
%    Ec2*ones(1,b/L*N),Ec1*ones(1,ceil(a*N/L))];
% Ev=[Ev1*ones(1,ceil(a*N/L)),Ev2*ones(1,b/L*N),Ev1*ones(1,ceil(a*N/L)),...
%    Ev2*ones(1,b/L*N),Ev1*ones(1,ceil(a*N/L))];
Dn=@(E,me,Ec)sqrt(2*me.^3.*(E-Ec))/(pi^2*hbar^2);
Dp=@(E,mh,Ev)sqrt(2*mh.^3.*(Ev-E))/(pi^2*hbar^2);
fED=@(E,Ef)1./(1+exp((E-Ef)/(k*T)));

ro1=-integral(@(E)Dn(E,m1e,Ec1).*(fED(E,Ec1/2)),Ec1,150*k*T)+integral(@(E)Dp(E,m1h,Ev1).*(fED(E,Ec1/2)),-150*k*T,Ev1);
ro2=-integral(@(E)Dn(E,m2e,Ec2).*(fED(E,Ec2/2)),Ec2,150*k*T)+integral(@(E)Dp(E,m2h,Ev2).*(fED(E,Ec2/2)),-150*k*T,Ev2);

ro=e*[ro1*ones(1,ceil(a*N/L)),ro2*ones(1,b/L*N),ro1*ones(1,ceil(a*N/L)),...
    ro2*ones(1,b/L*N),ro1*ones(1,ceil(a*N/L))];
