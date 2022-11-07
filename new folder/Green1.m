
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
U0=1;                   %height of barrier, 1 eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*0.3)*m0;%eff. mass in AlGaAs(0.3), kg
Ef=0.14;                %Fermi energy, (at T=300)
T=300;                  %room temperature, K

L=11e-9;
Np=110;                        %amount of steps
x=linspace(0,L, Np);       %creating a 'x-axis'
dx=x(2)-x(1);                   %definig a primive step
koef=-hbar^2/(2*m0*(dx^2))/e;   %definig a coefficient for analytical solving


% E=eye(Np)*(30);
% E=E+diag(ones(1,Np-1)*(-16),-1);
% E=E+diag(ones(1,Np-1)*(-16),1);
% E=E+diag(ones(1,Np-2)*(1),-2);
% E=E+diag(ones(1,Np-2)*(1),2);
E=eye(Np)*(-2);
E=E+diag(ones(1,Np-1)*(1),-1);
E=E+diag(ones(1,Np-1)*(1),1);


E=E*koef;
E(:,1                    :floor((b)/L*Np))    =E(:,1                    :floor((b)/L*Np))/0.067;
E(:,floor((b)/L*Np)+1    :floor((2*b)/L*Np))  =E(:,floor((b)/L*Np)+1    :floor((2*b)/L*Np))/0.0919;
E(:,floor((2*b)/L*Np)+1  :floor((2*b+a)/L*Np))=E(:,floor((2*b)/L*Np)+1  :floor((2*b+a)/L*Np))/0.067;
E(:,floor((2*b+a)/L*Np)+1:floor((3*b+a)/L*Np))=E(:,floor((2*b+a)/L*Np)+1:floor((3*b+a)/L*Np))/0.0919;
E(:,floor((3*b+a)/L*Np)+1:floor((4*b+a)/L*Np))=E(:,floor((3*b+a)/L*Np)+1:floor((4*b+a)/L*Np))/0.067;
U=zeros(1,Np);
U(floor((b)/L*Np):floor((2*b)/L*Np)-1)=U0;
U(floor((2*b+a)/L*Np):floor((3*b+a)/L*Np)-1)=U0;
H=E+diag(U);

Umax=1.1;
Ener=linspace(0,Umax,5*Np);
D=zeros(1,5*Np);

for i=1:5*Np
En=Ener(i);
oneone=zeros(Np);   oneone(1,1)=1;
lastlast=zeros(Np); lastlast(Np,Np)=1;
S1=zeros(Np)+oneone*koef.*exp(1i*sqrt(2*m1*En*e)/hbar*dx);
S2=zeros(Np)+lastlast*koef.*exp(1i*sqrt(2*m1*En*e)/hbar*dx);

G1=1i*(S1-S1');
G2=1i*(S2-S2');

Gr=(eye(Np)*(En+1i*1e-12)-H-S1-S2);

D(i)=real(trace(G1/Gr*G2/Gr'));
end

figure ('Units','normalized','OuterPosition', [0 0 1 1])
subplot(1,3,2)
semilogx(D,Ener)
ylabel('$E,eV$','Interpreter','latex')
xlabel('$log(D)$','Interpreter','latex')
title('$Transmission$ $coefficient$ $(log)$','Interpreter','latex')
ylim([0, Umax])

subplot(1,3,3)
plot(D,Ener)
ylabel('$E,eV$','Interpreter','latex')
xlabel('$D$','Interpreter','latex')
title('$transmission$ $coefficient$','Interpreter','latex')
ylim([0, Umax])
xlim([0 1])

subplot(1,3,1)
plot(x*1e9,U)
hold on
plot(x*1e9,(U+67/30)*0.03)
ylim([0, Umax])
xlim([x(1), x(end)]*1e9)
xlabel('$x,nm$','Interpreter','latex')
ylabel('$U, eV$','Interpreter','latex')
title('$Structure$','Interpreter','latex')






