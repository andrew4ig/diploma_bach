%preparing workspace
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
U0=1;                 %height of barrier, 1 eV
m1=0.067*m0;             %eff. mass in GaAs, kg
m2=(0.067+0.083*0.3)*m0;%eff. mass in AlGaAs(0.3), kg
Ef=0.14*e;              %Fermi energy, (at T=300)
% T=300;                  %room temperature, K
% v=1e6;                  %frequency, 1Mhz
% Ez=0.05*e;              %connection energy
% V=@(t)sin(2*pi*v*t);    %periodic voltage, V

L=11e-9;
Np=110;                        %amount of steps
x=linspace(-L/2,L/2, Np);       %creating a 'x-axis'
dx=x(2)-x(1);                   %definig a primive step
koef=-hbar^2/(2*m0*(dx^2))/(e);   %definig a coefficient for numeric solving

% 
% E=eye(Np)*(30);
% E=E+diag(ones(1,Np-1)*(-16),-1);
% E=E+diag(ones(1,Np-1)*(-16),1);
% E=E+diag(ones(1,Np-2)*(1),-2);
% E=E+diag(ones(1,Np-2)*(1),2);
E=eye(Np)*(-2);
E=E+diag(ones(1,Np-1)*(1),-1);
E=E+diag(ones(1,Np-1)*(1),1);
E=E*koef;
E(:,1                    :floor((b)/L*Np)-1)    =E(:,1                    :floor((b)/L*Np)-1)/0.067;
E(:,floor((b)/L*Np)      :floor((2*b)/L*Np)-1)  =E(:,floor((b)/L*Np)      :floor((2*b)/L*Np)-1)/0.0919;
E(:,floor((2*b)/L*Np)    :floor((2*b+a)/L*Np))  =E(:,floor((2*b)/L*Np)    :floor((2*b+a)/L*Np))/0.067;
E(:,floor((2*b+a)/L*Np+1):floor((3*b+a)/L*Np))  =E(:,floor((2*b+a)/L*Np+1):floor((3*b+a)/L*Np))/0.0919;
E(:,floor((3*b+a)/L*Np+1):floor((4*b+a)/L*Np))  =E(:,floor((3*b+a)/L*Np+1):floor((4*b+a)/L*Np))/0.067;
U=zeros(1,Np);
U(:,floor((b)/L*Np):floor((2*b)/L*Np)-1)=U0;
U(:,floor((2*b+a)/L*Np+1):floor((3*b+a)/L*Np))=U0;
H=E+diag(U);

oneone=zeros(Np);   oneone(1,1)=1;
lastlast=zeros(Np); lastlast(Np,Np)=1;
S1=@(En)zeros(Np)-oneone*koef.*exp(1i*sqrt(2*m1*En*e)/hbar*dx);
S2=@(En)zeros(Np)-lastlast*koef.*exp(1i*sqrt(2*m1*En*e)/hbar*dx);
% S1=@(E)zeros(Np)+oneone.*E/(2*pi*m1);
% S2=@(E)zeros(Np)+lastlast.*E/(2*pi*m1);

G1=@(En)1i*(S1(En)-S1(En)');
G2=@(En)1i*(S2(En)-S2(En)');

Gr=@(En)inv(eye(Np).*(En+1i*1e-12)-H-S1(En)-S2(En));

D=@(En)(trace(G1(En)*Gr(En)*G2(En)*(Gr(En)')));

Ener=linspace(0,1*e,Np);
% T=D(Ener);
% 
T=zeros(1,Np);
for i=1:Np
   T(i)=D(Ener(i));
end
plot(Ener/e, T)




