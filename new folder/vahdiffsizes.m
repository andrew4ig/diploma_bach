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
a=3e-9;                 %width of well
b=2e-9;                 %width of barrier
U0=1;                   %height of barrier, 1 eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*0)*m0*0+m1;%eff. mass in AlGaAs(0.3), kg
n=1e21;
Nen=300;
V=linspace(0,1,Nen);
V3=0;

%defining structure
j=2;                    %amount of barriers
dx=2e-10;               %grid step
L=j*b+(j+1)*a;          %total length
Np=floor(L/dx);         %amount of grid cells
koef=hbar^2/(2*m1*(dx^2))/e;
x=(0:Np-1)*dx;          %x-vector
One=linspace(0,1,Np);   %empty vector
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.2; %chem potential
% mu=mu*(1-pi^2/12*(k*300/e/mu)^2)

%correting hamiltonian corresponding to the heterstructure
E=eye(Np)*(2);
E=E+diag(ones(1,Np-1)*(-1),-1);
E=E+diag(ones(1,Np-1)*(-1),1);
E=E*koef;
for t=0:(2*j-1)
    if(mod(t,2)==0)
        E(:,round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))=E(:,round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))/(m1/m2);
    end
end
U=zeros(1,Np);
for t=0:(2*j-1)
    if(mod(t,2)==0)
        U(round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))=U0;
    end
end
H=E+diag(U);

Umax=1.1;
Ener=linspace(0,Umax,Nen);

D=zeros(1,Nen);
J=zeros(1,length(V));
JJ=zeros(1,length(V));
for klo=1:length(V)
    U1=-V(klo)*One;
    Ham=H+diag(U1);
    S1=zeros(Np);S2=S1;
    for i=1:Nen
        En=Ener(i);
        S1(1,1)=-koef*exp(1i*sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
        S2(Np,Np)=-koef*exp(1i*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);

        G1=1i*(S1-S1');
        G2=1i*(S2-S2');

        Gr=(eye(Np)*(En)-Ham-S1-S2);
        D(i)=real(trace(G1/Gr*G2/Gr'));
    end
   
    %transmission coefficient w/o votalge
    if(klo==1)
        Emax=Ener(islocalmax(D));
        D1=D;
        if (isempty(Emax))
            Emax=0;
            D1=0;
        end
    end
    
    %fermi functions for source and drain
    f1=1./(1+exp(-(mu-Ener)./(k*T/e)));
    f2=1./(1+exp(-(mu-V(klo)-Ener)./(k*T/e)));
    
    dE=Ener(2)-Ener(1); %energy step
     
    %J by Green
    J(klo)=e/(2*pi*hbar)*trapz(D.*(f1-f2))*dE*e;

    %J by Tsu
    S=@(Ez)log((1+exp((mu-Ez)/(k*T/e)))./(1+exp((mu-Ez-V(klo))/(k*T/e))));
    S=S(Ener);
    JJ(klo)=m1*e*k*T/(2*pi^2*hbar^3/e)*trapz(S.*D)*dE*e;
end

%microcurrent
J=J*1e6;
JJ=JJ*1e6;

f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
subplot(1,3,2)
plot(D1,Ener)
ylabel('E,????')
xlabel('D')
title('?????????????????????? ??????????????????????')
ylim([0, Umax])
xlim([0 1])
hold on
% y1=yline(mean(Emax(1)), '-',[' ',num2str(round(mean(Emax(1))*1000)), ' ??????']);
% y1.LabelHorizontalAlignment = 'center';

subplot(1,3,1)
plot(x*1e9,U)
hold on
ylim([-V3, Umax-V3])
xlim([x(1), x(end)]*1e9)
xlabel('x,????')
ylabel('U, ????')
title('???????????????????????????? ?????????????? ????????????????')

% shift=JJ./J;
% shift=mean(shift(2:end))
% 1/shift
subplot(1,3,3)
hold off
plot(V,J)
hold on
plot(V,JJ,'--k','LineWidth',2)
ylabel('J, ????A')
xlabel('V, ????')
title('??????')
legend('green','Tsu')
% exportgraphics(f,['2vcc',num2str(V3),'.jpg'])
    
