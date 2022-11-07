%preparing workspace
clc
clear
% close all


%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
T=300;

% X=0.3;
% %initial paramets
% a=3e-9;                 %size of chanel, nm
% b=2e-9;                 %size of barrier, nm
U0=0.6749;
m1=0.041*m0;
m2=(0.067+0.083)*m0;

Nen=300;
V=linspace(0,1,Nen);
V3=0;

%defining structure
j=2;                    %amount of barriers
dx=1e-10;               %grid step
L=13.7e-9;
N1=floor(L/dx);
x=(0:N1-1)*dx;          %total length
koef=hbar^2/(2*m1*(dx^2))/e;


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

Umax=0.7;%1*X;
Ener=linspace(0,Umax,Nen);

D=zeros(1,Nen);
J=zeros(1,length(V));
JJ=zeros(1,length(V));

% v = VideoWriter('Green1eV.avi','Motion JPEG AVI');
% v.Quality = 95;
% v.FrameRate=5;
% open(v);


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

        Gr=(eye(Np)*(En)-Ham-S1-S2);
        D(i)=real(trace(G1/Gr*G2/Gr'));
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
    
    dE=Ener(2)-Ener(1); %energy step


    %J by Tsu
    S=@(Ez)log((1+exp((mu-Ez)/(k*T/e)))./(1+exp((mu-Ez-V(klo))/(k*T/e))));
    S=2*pi*m1*k*T/(hbar^2)*S(Ener);
    JJ(klo)=e/hbar*2/(2*pi)^3*trapz(S.*D)*dE*e;
end



%microcurrent
% J=J*1e6;
% JJ=JJ*1e6;
% 
f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

subplot(1,3,2)
plot(D1,Ener)
ylabel('E,эВ')
xlabel('D')
title('Коэффициент прохождения')
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
title('Энергетический профиль стуктуры')

subplot(1,3,3)
hold off
% plot(V,J*e*4)
hold on
plot(V,JJ*e*4,'k','LineWidth',2)
ylabel('J, мкА/м^2')
xlabel('V, эВ')
title('ВАХ')
grid on
legend('green','Tsu','Location','best')

