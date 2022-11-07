%preparing workspace
clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
T=300;

%initial paramets
a=2e-9;                 %size of chanel, 2nm
b=2e-9;                 %size of barrier, 3nm
U0=1;                   %height of barrier, 1 eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*0.3)*m0;%eff. mass in AlGaAs(0.3), kg
n=1e15;
% N=10;                   %amount of iterations
Nen=800;
V=linspace(0,15,100);
V3=0;

j=2;
dx=2e-10;
L=(j+1)*a+j*b;
Np=floor(L/dx);
koef=-hbar^2/(2*m1*(dx^2))/e;
x=(0:Np-1)*dx;
One=linspace(0,1,Np);
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.2;
% mu=mu*(1-pi^2/12*(k*300/e/mu)^2)

E=eye(Np)*(-2);
E=E+diag(ones(1,Np-1)*(1),-1);
E=E+diag(ones(1,Np-1)*(1),1);
E=E*koef;
for t=1:j
    E(:,floor(((2*t-1)*b)/L*Np):floor((2*t*b)/L*Np))=E(:,floor(((2*t-1)*b)/L*Np):floor((2*t*b)/L*Np))/(0.067/0.0919);
end
U=zeros(1,Np);
if(j==3)
    U=-V3*ones(1,Np)+[V3*ones(1,10),(1-linspace(0,1,39))*V3,zeros(1,21)];
end
for t=1:j
    U(floor(((2*t-1)*b)/L*Np):floor((2*t*b)/L*Np))=U(floor(((2*t-1)*b)/L*Np):floor((2*t*b)/L*Np))+U0;
end
H=E+diag(U);

Umax=1.5;
Ener=linspace(0,Umax,Nen);

D=zeros(1,Nen);

J=zeros(1,length(V));
for klo=1:length(V)
    klo
    U1=U-V(klo)*One;
    Ham=H+diag(U1);
    S1=zeros(Np);S2=S1;
    for i=1:Nen
        En=Ener(i);
        S1(1,1)=koef*exp(1i*sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
        S2(Np,Np)=koef*exp(1i*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);

        G1=1i*(S1-S1');
        G2=1i*(S2-S2');

        Gr=(eye(Np)*(En)-Ham-S1-S2);
        D(i)=real(trace(G1/Gr*G2/Gr'));
    end
   
    if(klo==1)
        Emax=Ener(islocalmax(D));
        D1=D;
    end
    
    f1=1./(1+exp(-(mu-Ener)./(k*T/e)));
    f2=1./(1+exp(-(mu-V(klo)-Ener)./(k*T/e)));
    
    dE=Ener(2)-Ener(1);
    I=sum(D.*(f1-f2))*dE;
    J(klo)=e/(2*pi*hbar)*I;
end

J=J*e*1e6;

f=figure ('Units','normalized','OuterPosition', [0 0 0.4 0.4]);
subplot(1,3,2)
plot(D1,Ener)
ylabel('E,эВ')
xlabel('D')
title('Коэффициент прохождения')
ylim([0, Umax])
xlim([0 1])
hold on
y1=yline(mean(Emax(1)), '-',[' ',num2str(round(mean(Emax(1))*1000)), ' мэВ']);
y1.LabelHorizontalAlignment = 'center';


subplot(1,3,1)
plot(x*1e9,U)
hold on
ylim([-V3, Umax-V3])
xlim([x(1), x(end)]*1e9)
xlabel('x,нм')
ylabel('U, эВ')
title('Энергетический профиль стуктуры')

subplot(1,3,3)
plot(V,J)
hold on
%     
%     xlim([x(1), x(end)]*1e9)
ylabel('J, мкA')
xlabel('V, эВ')
title('ВАХ')

exportgraphics(f,['2vcc',num2str(V3),'.jpg'])
    
