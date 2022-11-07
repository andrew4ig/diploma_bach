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

X=0.3;
%initial paramets
a=3e-9;                 %size of chanel, nm
b=2e-9;                 %size of barrier, nm
U0=X*0.74;            %height of barrier,  eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg

Nen=500;
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
n=1e21;
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.3; %chem potential
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

% Umax=1.1;
Umax=1*X;
Ener=linspace(0,Umax,Nen);

D=zeros(1,Nen);
J=zeros(1,length(V));
JJ=zeros(1,length(V));

for klo=1:length(V)

    %correcting hamiltonian
    U1=-V(klo)*One;
    Ham=H+diag(U1);

    %Green's transmition coef-t
    S1=zeros(Np);S2=S1;
    for i=1:Nen
        En=Ener(i);
        S1(1,1)=-koef*exp(1i*sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
        S2(Np,Np)=-koef*exp(1i*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);

        G1=1i*(S1-S1');
        G2=1i*(S2-S2');

        Gr=(eye(Np)*(En)-Ham-S1-S2); %инициализация функции Грина
        D(i)=real(trace(G1/Gr*G2/Gr')); %расчет прозрачности
    end
   
    %transmission coefficient w/o votalge
    if(klo==1)
        Emax=Ener(islocalmax(D));
        D1=D;
        if (isempty(D1))
            Emax=0;
            D1=0;
        end
    end
    
    %fermi functions for source and drain
    f1=@(Ex,Ey,Ez)1./(1+exp((Ex+Ey+Ez-mu)./(k*T/e)));
    f2=@(Ex,Ey,Ez)1./(1+exp((Ex+Ey+Ez-mu+V(klo))./(k*T/e)));
    
    dE=Ener(2)-Ener(1); %energy step
    
     %attempt to integrate numerically
%     Uup=10; Steps=105;
%     ex=linspace(0,Uup,Steps);
%     ey=linspace(0,Uup,Steps);
%     [EX,EY,EZ]=meshgrid(ex,ey,Ener);
% 
%     f=f1(EX,EY,EZ)-f2(EX,EY,EZ);
%     
%    %F=  squeeze(trapz(ey,trapz(ex,f,1),2));

    %inegrating functionally
    q1=@(x)integral(@(y)(f1(x,y,Ener)-f2(x,y,Ener)),0,Inf,'RelTol',1e-8,'AbsTol',1e-13,'ArrayValued',true);
    F=integral(q1,0,Inf,'RelTol',1e-8,'AbsTol',1e-13,'ArrayValued',true)'*10;
  
%     %fermi functions for source and drain
%     f1=1./(1+exp(-(mu-Ener)./(k*T/e)));
%     f2=1./(1+exp(-(mu+V(klo)-Ener)./(k*T/e)));
% 

    %J by Green
    J(klo)=e/hbar*2/(2*pi)^3*trapz(D.*(F'))*dE;

    %J by Tsu
    S=@(Ez)log((1+exp((mu-Ez)/(k*T/e)))./(1+exp((mu-Ez-V(klo))/(k*T/e))));
    S=2*pi*m1*k*T/hbar^2*S(Ener)*e;
    JJ(klo)=e/hbar*2/(2*pi)^3*trapz(S.*D)*dE;
   
    if(klo==1)
        h=figure('Units','normalized','OuterPosition', [0 0 0.5 1]);
        subplot(3,1,1)
        plot(Ener, log10(S),'--k','LineWidth',2);
        hold on;
        plot(Ener, log10(F))
        legend ('S','int')
        grid on
        xlabel('Ez, эВ')
        ylabel('Ф-ция снаб-ия, lg')

        subplot(3,1,2)
        plot(Ener, (S),'--k','LineWidth',2);
        hold on;
        plot(Ener, (F))
%         plot(Ener, 1./(S').*(F));
%         plot(Ener, (S')-(F));
        legend ('S','int')%,'/','-')
        grid on
        xlabel('Ez, эВ')
        ylabel('Ф-ция снаб-ия, б/р')
        yyaxis right;
        plot(Ener, D,'-.m','HandleVisibility','off')

        subplot(3,1,3)
        plot(Ener, (S.*D),'--k','LineWidth',2);
        hold on;
        plot(Ener, (F'.*D))
        xlabel('Ez, эВ')
        ylabel('S*D')
        grid on
    end
end

%%
f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.4 0.45]);
subplot(1,3,2)
plot(D1,Ener)
ylabel('E,эВ')
xlabel('D')
title({'Коэффициент', 'прохождения'})
set(gca,'FontSize',11,'fontWeight','bold')
ylim([0, U0+0.1])
xlim([0 1])
grid on
hold on

subplot(1,3,1)
plot(x*1e9,U)
hold on
ylim([-V3, U0+0.1-V3])
xlim([x(1), x(end)]*1e9)
xlabel('x,нм')
ylabel('U, эВ')
grid on
title({'Энергетический профиль' ,'стуктуры'})
set(gca,'FontSize',11,'fontWeight','bold')

subplot(1,3,3)
hold off
plot(V,J)
hold on
plot(V,JJ,'--k','LineWidth',2)
ylabel('J, А/м^2')
xlabel('V, эВ')
title('ВАХ')
grid on
legend('f_1-f_2','Ф. снабжения','Location','best')
set(gca,'FontSize',11,'fontWeight','bold')

%  exportgraphics(h,[num2str(number),'supply.jpg']) 
%  exportgraphics(f,[num2str(number),'vcc.jpg']) 
% number=number+1