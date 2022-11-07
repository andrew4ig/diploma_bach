%preparing workspace
clc
clear
% close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;

X=0.3;
%initial paramets
a=3e-9;                 %size of chanel, nm ?????????????
b=2e-9;                 %size of barrier, nm ?????????????
U0=e*X*0.74;            %height of barrier,  eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg

g=@(E,m,U)sqrt(2*m*(E-U))/hbar;
g1=@(E)g(E,m1,0); g2=@(E)g(E,m2,U0);
g={@(E)g1(E);@(E)g2(E);@(E)g1(E);@(E)g2(E);@(E)g1(E)};
z=[0, b, b+a, 2*b+a]+b;
m=[m1,m2,m1,m2,m1];

T=@(k,E)[0.5*(1+g{k}(E)/g{k+1}(E)*m(k+1)/m(k))*exp(-1i*(g{k+1}(E)-g{k}(E))*z(k)),   0.5*(1-g{k}(E)/g{k+1}(E)*m(k+1)/m(k))*exp(-1i*(g{k+1}(E)+g{k}(E))*z(k));
       0.5*(1-g{k}(E)/g{k+1}(E)*m(k+1)/m(k))*exp(1i*(g{k+1}(E)+g{k}(E))*z(k)),      0.5*(1+g{k}(E)/g{k+1}(E)*m(k+1)/m(k))*exp(1i*(g{k+1}(E)-g{k}(E))*z(k))];
T0=@(E)T(4,E)*T(3,E)*T(2,E)*T(1,E);
Tline=@(E)reshape(T0(E),1,[]);
Umax=0.7*e;%e*X;
En=linspace(0.01*e, Umax,500);
D=zeros(1,500);
for i=1:500
    T1=Tline(En(i));
    D(i)=abs(g{5}(En(i)))/abs(g{1}(En(i)))*m(1)/m(5)*abs((T1(4)*T1(1)-T1(2)*T1(3))/T1(4))^2;
end

% figure('Units','normalized','OuterPosition',[0 0 1 1])

% subplot(1,3,2)
hold on
plot(D,En/e,'--k','LineWidth',1.5);
Emax=En(islocalmax(D)); 
% yline(Emax(1)/e, '-.',['first peak ', num2str(round(Emax(1)/e*1000)), ' meV'],'HandleVisibility','off');
% yline(Emax(2)/e, '-.',['second peak ', num2str(round(Emax(2)/e*1000)), ' meV']);
ylabel('E,эВ')
xlabel('D')
title('Коэффициент прохождения')

% subplot(1,3,3)
% semilogx(D,En/e);
% Emax=En(islocalmax(D)); 
% yline(Emax(1)/e, '-.',['first peak ', num2str(round(Emax(1)/e*1000)), ' meV']);
% % yline(Emax(2)/e, '-.',['second peak ', num2str(round(Emax(2)/e*1000)), ' meV']);
% title('$Transmission$ $coefficient$ $log$','Interpreter','latex')
% xlabel('$D$','Interpreter','latex')
% ylabel('$E,  ev$','Interpreter','latex')
% 
% 
% subplot(1,3,1)
% x=[0,a,a,a+b,a+b,2*a+b,2*a+b,2*b+2*a,2*b+2*a,3*a+2*b]*1e9;
% y=[0,0, U0, U0, 0, 0, U0, U0, 0, 0];
% plot(x,y/e)
% ylim([0, Umax/e])
% xlim([0 x(end)])
% title('$Structure$','Interpreter','latex')
% xlabel('$x,nm$','Interpreter','latex')
% ylabel('$E,  ev$','Interpreter','latex')
% yline(Emax(1)/e, '-.',['first peak ', num2str(round(Emax(1)/e*1000)), ' meV']);
% % yline(Emax(2)/e, '-.',['second peak ', num2str(round(Emax(2)/e*1000)), ' meV']);
