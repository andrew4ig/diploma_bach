%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
T=300;

%initial paramets
a=3e-9;                 %size of chanel, nm
b=2e-9;                 %size of barrier, nm
X=0.3;                  %replacement proportion
U0=X*0.74;              %height of barrier,  eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg

%defining structure
j=2;                    %amount of barriers
dx=2e-10;               %grid step
L=j*b+(j+1)*a;          %total length
Np=floor(L/dx);         %amount of grid cells
x=(0:Np-1)*dx;          %x-vector
koef=hbar^2/(2*m1*(dx^2))/e;


%correting hamiltonian corresponding to the heterostructure
E=eye(Np)*(2);
E=E+diag(ones(1,Np-1)*(-1),-1);
E=E+diag(ones(1,Np-1)*(-1),1);
E=E*koef;
for t=0:(2*j-1)
    if(mod(t,2)==0)
        E(:,round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))=...
        E(:,round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))/(m1/m2);
    end
end
U=zeros(1,Np);
for t=0:(2*j-1)
    if(mod(t,2)==0)
        U(round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))=U0;
    end
end
H=E+diag(U);

Nen=500;                %amount of voltage steps
Vmax=1;                 %maximun voltage
V=linspace(0,1,Nen);    %voltage vector
Umax=1*X;               %maximum energy
Ener=linspace(0,Umax,Nen);%energy vector
One=linspace(0,1,Np);   %empty vector
n=3e21; dE=Ener(2)-Ener(1);
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.3; %chem potential

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

    S=@(Ez)2*pi*m1*k*T/hbar^2*log((1+exp((mu-Ez)/(k*T/e)))./(1+exp((mu-Ez-V(klo))/(k*T/e))));
    
    %current
    J(klo)=e/hbar*2/(2*pi)^3*trapz(D.*(S(Ener)))*(dE*e);

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
ylabel('J, А/м^2')
xlabel('V, эВ')
title('ВАХ')
grid on
set(gca,'FontSize',11,'fontWeight','bold')

rectangle('Position',[0.3 2.1e11 0.1 0.25e11])
axes('Position',[0.775 0.17 0.12 0.3])
box on
indexs = V>0.3 & V<0.4;
hold on
plot(V(indexs),JJ(indexs))
plot(V(indexs),J(indexs),'k--','LineWidth',2)
% set(gca, 'YTick', []);
% set(gca, 'XTick', []);
grid minor

% JJ – ток по Цу-Есаки
% J - ток по Грину
max((JJ-J)./J*100)
