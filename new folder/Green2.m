%preparing workspace
clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
n=1e15;
T=300;

%initial paramets
a=2e-9;                 %size of chanel, 3nm
b=2e-9;                 %size of barrier, 2nm
U0=1;                   %height of barrier, 1 eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*0.3)*m0;%eff. mass in AlGaAs(0.3), kg
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.2;
% mu=mu*(1-pi^2/12*(k*300/e/mu)^2)

L=10e-9;
Np=100;                        %amount of steps
x=linspace(0,L, Np);            %creating a 'x-axis'
dx=x(2)-x(1);                   %definig a primive step
koef=-hbar^2/(2*m1*(dx^2))/e;   %definig a coefficient for analytical solving

E=eye(Np)*(-2);
E=E+diag(ones(1,Np-1)*(1),-1);
E=E+diag(ones(1,Np-1)*(1),1);
E(:,floor((b)/L*Np)+1    :floor((2*b)/L*Np))  = E(:,floor((b)/L*Np)+1    :floor((2*b)/L*Np))/(0.067/0.0919);
E(:,floor((2*b+a)/L*Np)+1:floor((3*b+a)/L*Np))= E(:,floor((2*b+a)/L*Np)+1:floor((3*b+a)/L*Np))/(0.067/0.0919);
E=E*koef;

U=zeros(1,Np);
U(floor((b)/L*Np):floor((2*b)/L*Np))=U0;
U(floor((2*b+a)/L*Np):floor((3*b+a)/L*Np))=U0;
H=E+diag(U);

Multipl=2;
Umax=1.5;
V=linspace(0,Umax);
One=linspace(0,1,Np);

Ener=linspace(0,Umax,Multipl*Np);
D=zeros(1,Multipl*Np);
S1=zeros(Np);S2=S1;
J=zeros(1,Np);
for j=1:length(V)
    j
    U1=U-V(j)*One;
    Ham=H+diag(U1);
    for i=1:Multipl*Np
        En=Ener(i);
        S1(1,1)=koef*exp(1i*sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
        S2(Np,Np)=koef*exp(1i*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);

        G1=1i*(S1-S1');
        G2=1i*(S2-S2');

        Gr=(eye(Np)*(En)-Ham-S1-S2);
        D(i)=real(trace(G1/Gr*G2/Gr'));
    end
    if(j==1)
        D1=D;
        Emax1=Ener(islocalmax(D)); 
    end 
  
    f1=1./(1+exp(-(mu-Ener)./(k*T/e)));
    f2=1./(1+exp(-(mu-V(j)-Ener)./(k*T/e)));
    
    dE=Ener(2)-Ener(1);
    I=sum(D.*(f1-f2))*dE;
    J(j)=e/(2*pi*hbar)*I;
end

J=J*e;
figure('Units','normalized','OuterPosition',[0 0 1 1])

subplot(1,3,1)

y=(U);
plot(x,y)
ylim([0, 1.5])
xlim([0 x(end)])
title('Профиль структуры')
xlabel('x, нм')
ylabel('E, эВ')
yline(Emax1(1), '-.',[num2str(round(Emax1(1)*1000)), ' мВ']);

subplot(1,3,2)
plot(D1,Ener);
yline(Emax1(1), '-.',[num2str(round(Emax1(1)*1000)), ' мВ']);
title('Коэффициент прохождения')
xlabel('D')
ylabel('E, эВ')

subplot(1,3,3)
plot(V,J*1e6)
hold on
ylabel('J, мкA')
xlabel('V, эВ')
title('ВАХ')





