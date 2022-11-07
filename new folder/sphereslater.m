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
T=300;                  %room temperature, K
a=3e-9;                 %size of chanel, 3nm
b=2e-9;                 %size of barrier, 2nm
U0=1*e;                 %height of barrier, 1 eV
m=0.067*m0;             %eff. mass in GaAs, kg
mb=(0.067+0.083*0.5)*m0;%eff. mass in AlGaAs(0.5), kg
Ef=0.14*e;              %Fermi energy, eV (at T=300)(0,14)
Eg=1.42*e;              %gap energy, eV
v=1e6;                  %frequency, 1Mhz
Ez=0.1*e;             %connection energy
V=@(t)1*sin(2*pi*v*t);    %periodic voltage, V

eps0=0.244*e;

% eps=@(nw,t)(eps0-e*V(t)/2+Ez*nw+abs(eps0-e*V(t)/2+Ez*nw))/2;
eps=@(nw,t)(eps0-e*V(t)/2+Ez*nw);

Eks=@(nw,t)eps(nw,t);
Ekd=@(nw,t)eps(nw,t)+e*V(t);
Ekw=@(nw,t)eps0;

wsw=@(nw,t)sqrt(2*Eks(nw,t)/m)/b;
wdw=@(nw,t)sqrt(2*Ekd(nw,t)/m)/b;
wws=@(nw,t)sqrt(2*Ekw(nw,t)/m)/b;
wwd=@(nw,t)wws(nw,t);

f2d=@(EF,E,E0) log(1+exp(-(E-EF)/(k*T))).*heaviside(E-E0);

% E=linspace(0,2)*e;
% figure('Units','normalized','OuterPosition',[.2 .2 .6 .6])
% plot(E/e,f2d(Ef,E,0))
% xline(Eg/e,'-','E_c')
% xline(Ef/e,'--','E_f')
% xline(0/e,'--','E_v')
% xline(Ef/e-1,'--','E_f-eV')
% xline(Ef/e+1,'--','E_f+eV')
% xline(eps0/e, '--','\epsilon_0')
% grid on
% ylabel('$E,eV$','Interpreter','latex')
% xlabel('$f_{2d}$','Interpreter','latex')
% title('$2D$ $distribution$','Interpreter','latex')
% close all

Efs=Ef;
Efd=@(t)Ef-e*(V(t));

ns=@(nw,t)f2d(Efs,eps(nw,t),0);       %
nd=@(nw,t)f2d(Efd(t),eps(nw,t),-e*(V(t)));    %


jsw=@(nw,t)wsw(nw,t).*ns(nw,t).*(1-nw);
jdw=@(nw,t)wdw(nw,t).*nd(nw,t).*(1-nw);
jws=@(nw,t)wws(nw,t).*nw;
jwd=@(nw,t)wwd(nw,t).*nw;

eq = @(nw,t)jsw(nw,t)+jdw(nw,t)-jws(nw,t)-jwd(nw,t);
Nw0=0;
[t,Nw]=ode15s(@(t,nw)eq(nw,t),[0 1/v],Nw0,odeset('Mass',10,'RelTol',1e-10));

J=jsw(Nw,t)-jdw(Nw,t)-jws(Nw,t)+jwd(Nw,t);
V=V(t);

% figure('Units','normalized','OuterPosition',[0 0 1 1])
x=[0,b,b,2*b,2*b,2*b+a,2*b+a,3*b+a,3*b+a,4*b+a]*1e9;
y=[0 0 1 1 0 0 1 1 0 0];
Js=zeros(1,4);temp=1;
k=0;
for i=[75 107 126 299] 
    f1=figure('Units','normalized','OuterPosition',[0.1 0.1 0.25 0.25]);
    subplot(1,2,1)
    y=[0 0 1 1-b/(2*b+a)*V(i) 0-b/(2*b+a)*V(i) 0-(a+b)/(2*b+a)*V(i) 1-(a+b)/(2*b+a)*V(i) 1-V(i) 0-V(i) 0-V(i)];
    plot(x,y,'LineWidth',1.5)                                   %structure
    hold on
    grid on
    k=i;
    
    plot([2*b,2*b+a]*1e9,[1 1]*eps(Nw(k),t(k))/e,'LineWidth',1.5)   %well energy value
    text(1,1.3, sprintf('$V=%2.2fB$',V(k)),'Interpreter','latex','FontSize',12)

    J_1=jsw(Nw(k),t(k))-jdw(Nw(k),t(k))-jws(Nw(k),t(k))+jwd(Nw(k),t(k));
    if(k==75||k==107||k==126||k==299)
        Js(temp)=J_1;
        temp=temp+1;
    end

    Uf=linspace(min(-e*(V(k)),0), e);
    plot(-f2d(Efs,Uf,0)/4,Uf/e,'--','LineWidth',1.5)
    plot((4*b+a)*1e9+f2d(Efd(t(k)),Uf,-e*(V(k)))/4,Uf/e,'--','LineWidth',1.5)
    
    ylim([-1, 2])
    xlim([-4.5, 15.5])
    set(gca,'xtick',[])
    
    subplot(1,2,2)
    m=linspace(-1,1);
    n=@(t)sqrt(1-t.^2);
    M=[m,fliplr(m)]; N=[n(m),(-n(m))];
    plot(M,N);
    hold on
    
    drawArrow = @(x,y,varargin) quiver( 0,0,x,y,0, varargin{:} )  ;
    plot([0 0],[0 -1],'linewidth',2,'color','k')
    text(0.2,-0.5,'k_ф')
    k=sqrt((eps0-e*V(i)/2.1)/Ef);
    if (isreal(k))
        text(k/2,0.2,'k','FontWeight','bold')
        drawArrow(k,0,'linewidth',3,'color','k')
    else
        text(-0.5,0.50,'Im(k)<>0')
    end
    axis equal
    set(gca,'xtick',[],'ytick',[])
    exportgraphics(f1,['sphere_vcc',num2str(round(V(i)*100)/100),'.jpg'])
end

f2=figure('Units','normalized','OuterPosition',[0.7 0.5 0.3 0.5]);
plot(V,J*e,'LineWidth',1.5)
grid on
ylabel('J,{A / м^2*с}','FontSize',14)
xlabel('Напряжение,В','FontSize',14)
title('ВАХ','FontSize',16)
hold on
plot([V(75) V(113) V(133) V(299)],Js*e,'ro')
xline(V(75),'-.','V=0 В','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xline(V(113),'-.','V=0.28 В','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xline(V(133),'-.','V=0.47 В','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xline(V(299),'-.','V=0.54 В','LabelVerticalAlignment','top');
xlim([-0.1 1])
ylim([-0.1e-4 1.5e-4])
exportgraphics(f2,['sphere_vcc.jpg'])
