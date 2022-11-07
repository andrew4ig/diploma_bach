clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;

%initial paramets
a=3e-9;                 %size of chanel, 3nm
b=2e-9;                 %size of barrier, 2nm
U0=1*e;                 %height of barrieer, 1ev

%GaAs structure
eps1=12.90;
m1e=0.067*m0;
m1h=0.082*m0;
Eg1=1.42*e;
Ev1=0;
Ec1=Ev1+Eg1;
n1=1.1e13;

%GaAsAl0.05 structure
eps2=12.90-2.84*0.5;
m2e=(0.063+0.083*0.5)*m0;
m2h=(0.082+0.068*0.5)*m0;
Eg2=(1.9+0.125*0.5+0.143*0.5^2)*e;
Ev2=0;
Ec2=Ev2+Eg2;
n2=1.1e8;

T=300;                  %room temperature, K

L=3*a+2*b;  %structure length
N=L*1e9*10; %divideble number of steps
%x=linspace(-L/2,L/2,N);

%m_holes along x axis
mh=[m1h*ones(1,ceil(a*N/L)),m2h*ones(1,b/L*N),m1h*ones(1,ceil(a*N/L)),...
   m2h*ones(1,b/L*N),m1h*ones(1,ceil(a/L*N))];
%m_elec along x axis
me=[m1e*ones(1,ceil(a*N/L)),m2e*ones(1,b/L*N),m1e*ones(1,ceil(a*N/L)),...
   m2e*ones(1,b/L*N),m1e*ones(1,ceil(a*N/L))];

%energyes along x axis
Ec0=[Ec1*ones(1,ceil(a*N/L)),Ec2*ones(1,b/L*N),Ec1*ones(1,ceil(a*N/L)),...
   Ec2*ones(1,b/L*N),Ec1*ones(1,ceil(a*N/L))];
Ev0=[Ev1*ones(1,ceil(a*N/L)),Ev2*ones(1,b/L*N),Ev1*ones(1,ceil(a*N/L)),...
   Ev2*ones(1,b/L*N),Ev1*ones(1,ceil(a*N/L))];

%initial concetration along x axis
n=[n1*ones(1,ceil(a*N/L)),n2*ones(1,b/L*N),n1*ones(1,ceil(a*N/L)),...
   n2*ones(1,b/L*N),n1*ones(1,ceil(a/L*N))];
%eps along x axis
eps=[eps1*ones(1,ceil(a*N/L)),eps2*ones(1,b/L*N),eps1*ones(1,ceil(a*N/L)),...
   eps2*ones(1,b/L*N),eps1*ones(1,ceil(a/L*N))];

%supposin initial potential is zero
phi=0;

%correction by potential
Ec=Ec0-phi*e;
Ev=Ev0-phi*e;

%solving poission eq
phi=cumtrapz(cumtrapz(e*n./(eps0*eps)));

%solving Schr eq
m=m1e; mb=m2e; E0=0;
gm1=@(E)sqrt(2*m*(E-E0))/hbar;
gm2=@(E)sqrt(2*mb*(E-E0-U0))/hbar;
gm3=@(E)sqrt(2*m*(E-E0))/hbar;
gm4=@(E)sqrt(2*mb*(E-E0-U0))/hbar;
gm5=@(E)sqrt(2*m*(E-E0))/hbar;
A1=1;

MKoef=@(E)[ 1,         -1,                         -1,                                 0,                              0,                                  0,                                      0,                                      0;
            -gm1(E)/m, -gm2(E)/mb,                 gm2(E)/mb,                          0,                              0,                                  0,                                      0,                                      0;
            0,         exp(1i*gm2(E)*b),           exp(-1i*gm2(E)*b),                  -exp(1i*gm3(E)*b),              -exp(-1i*gm3(E)*b),                 0,                                      0,                                      0;
            0,         gm2(E)/mb.*exp(1i*gm2(E)*b), -gm2(E)/mb.*exp(-1i*gm2(E)*b),     -gm3(E)/m.*(1i*gm3(E)*b),       gm3(E)/m.*exp(-1i*gm3(E)*b),        0,                                      0,                                      0;
            0,         0,                          0,                                  exp(1i*gm3(E)*(a+b)),           exp(-1i*gm3(E)*(a+b)),              -exp(1i*gm4(E)*(a+b)),                  -exp(-1i*gm4(E)*(a+b)),                 0;
            0,         0,                          0,                                  gm3(E)/m.*exp(1i*gm3(E)*(a+b)), -gm3(E)/m.*exp(-1i*gm3(E)*(a+b)),   -gm4(E)/mb.*exp(1i*gm4(E)*(a+b)),       gm4(E)/mb.*exp(-1i*gm4(E)*(a+b)),       0;
            0,         0,                          0,                                  0,                              0,                                  exp(1i*gm4(E)*(a+2*b)),                 exp(-1i*gm4(E)*(a+2*b)),                exp(1i*gm5(E)*(a+2*b));
            0,         0,                          0,                                  0,                              0,                                  gm4(E)/mb.*exp(1i*gm4(E)*(a+2*b)),      -gm4(E)/mb.*exp(-1i*gm4(E)*(a+2*b)),   -gm5(E)/m.*exp(1i*gm5(E)*(a+2*b))];
MSvob=@(E)[-A1; -A1*gm1(E)/m; 0; 0; 0; 0; 0; 0];
Koef=@(E)MKoef(E)\MSvob(E);

Energyies=linspace(0,1.5,N)*e;
K=zeros(1,N);
for i=1:N
   temp=Koef(Energyies(i));
   K(i)=abs(temp(8)/A1)^2;
end

%resonant max
temp=islocalmax(K);
E=Energyies(temp);E=E(2);

%solving for a wave function
x=linspace(-a,L-a,N);
Ind=[A1,Koef(E)'];

x1=linspace(-a,0,ceil(a/L*N));
Psi1=Ind(1)*exp(1i*gm1(E)*x1)+Ind(2)*exp(-1i*gm1(E)*x1);
x2=linspace(0,b,ceil(b/L*N));
Psi2=Ind(3)*exp(1i*gm2(E)*x2)+Ind(4)*exp(-1i*gm2(E)*x2);
x3=linspace(b,b+a,ceil(a/L*N));
Psi3=Ind(5)*exp(1i*gm3(E)*x3)+Ind(6)*exp(-1i*gm3(E)*x3);
x4=linspace(b+a,a+2*b,ceil(b/L*N));
Psi4=Ind(7)*exp(1i*gm4(E)*x4)+Ind(8)*exp(-1i*gm4(E)*x4);
x5=linspace(a+2*b,2*b+2*a,ceil(a/L*N));
Psi5=Ind(9)*exp(1i*gm5(E)*x5);
Psi=[Psi1,Psi2,Psi3,Psi4,Psi5];
U=[zeros(1,ceil(a/L*N)),U0*ones(1,ceil(b/L*N)),zeros(1,ceil(a/L*N)),U0*ones(1,ceil(b/L*N)),zeros(1,ceil(a/L*N))];

% figure('Units','normalized','OuterPosition',[0 0 1 1])
% subplot(1,2,2)
% plot(K,Energyies/e)
% grid on
% subplot(1,2,1)
% plot(x,U/e)
% hold on
% plot(x,Psi/max(Psi)/5+E/e)
% grid on
% Psi=Psi/sqrt(Psi*Psi');
Psi=Psi*e*1000;

%using perturbation theory
dU=Psi*diag(phi)*Psi';

%enegry correction
Ea=E+dU;

%density of states and fermi distriburion
Dn=@(E,me,Ec)sqrt(2*me.^3.*(E-Ec))/(pi^2*hbar^2);
Dp=@(E,mh,Ev)sqrt(2*mh.^3.*(Ev-E))/(pi^2*hbar^2);
fED=@(E,Ef)1./(1+exp((E-Ef)/(k*T)));

% ro1=-integral(@(E)Dn(E,m1e,Ec1).*(fED(E,Ec1/2)),Ec1,150*k*T)+integral(@(E)Dp(E,m1h,Ev1).*(fED(E,Ec1/2)),-150*k*T,Ev1);
% ro2=-integral(@(E)Dn(E,m2e,Ec2).*(fED(E,Ec2/2)),Ec2,150*k*T)+integral(@(E)Dp(E,m2h,Ev2).*(fED(E,Ec2/2)),-150*k*T,Ev2);
% 
% ro=e*[ro1*ones(1,ceil(a*N/L)),ro2*ones(1,b/L*N),ro1*ones(1,ceil(a*N/L)),...
%     ro2*ones(1,b/L*N),ro1*ones(1,ceil(a*N/L))];
nn=integral(@(E)Dn(E,m1e,Ea/4).*(fED(E,Ea)),Ea,150*k*T);