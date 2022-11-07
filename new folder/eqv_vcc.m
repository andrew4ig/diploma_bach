%preparing workspace
clc
clear
close all

%constants
k=1.38e-23;
kT=k*300;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;

U0=0.6749;
m1=0.041*m0;
m2=(0.067+0.083)*m0;

dx1=1e-10;
L1=13.7e-9;
NE=1000;
N1=floor(L1/dx1);
x1=(0:N1-1)*dx1;

Vmax=0.3;
Volt=-Vmax:0.02:Vmax;
Umax=0.8*e;
Umax=max(Umax,U0*e);
En=linspace(0,Umax,NE); 

n=3e24;
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e;

g=@(E,m,U)sqrt(2*m*(E-U))/hbar;
T=@(m1,m2,g1,g2,z)[ 0.5*(1+g1/g2*m2/m1)*exp(-1i*(g2-g1)*z),   0.5*(1-g1/g2*m2/m1)*exp(-1i*(g2+g1)*z);
                    0.5*(1-g1/g2*m2/m1)*exp( 1i*(g2+g1)*z),   0.5*(1+g1/g2*m2/m1)*exp( 1i*(g2-g1)*z)];
J=zeros(1,length(Volt)); JJ=J;
D1=zeros(1,NE);

mm(1:N1)=m1;
mm(x1>=2e-9 & x1<4.5e-9)=m2;
mm(x1>=9.2e-9 & x1<11.7e-9)=m2;

v = VideoWriter('eqv.mp4','MPEG-4');
v.FrameRate = 2;
v.Quality = 100;
open(v);

for i=1:length(Volt)
    V=Volt(i);

    UU=zeros(1,N1);

    for ii=1:N1
        if(x1(ii)>=0&&x1(ii)<2e-9)
            UU(ii)=-V*ii/N1;
        elseif(x1(ii)>=2e-9&&x1(ii)<4.5e-9)
            UU(ii)=U0-V*ii/N1;
        elseif (x1(ii)>=4.5e-9&&x1(ii)<9.2e-9)
            UU(ii)=-V*ii/N1;
        elseif(x1(ii)>=9.2e-9&&x1(ii)<11.7e-9)
            UU(ii)=U0-V*ii/N1;
        elseif(x1(ii)>=11.7e-9&&x1(ii)<13.7e-9)
            UU(ii)=-V*ii/N1;
        end
    end

%      UU=(UU-V/N1/2)*e;
    UU=UU*e;


    T0=@(E)T(m2,m1,g(E,m2,UU(N1)),g(E,m1,-V*e),x1(N1));
    for ii=N1-1:-1:1
        T0=@(E)T0(E)*T(mm(ii),mm(ii+1),g(E,mm(ii),UU(ii)),g(E,mm(ii+1),UU(ii+1)),x1(ii));
    end
    T0=@(E)T0(E)*T(m1,m2,g(E,m1,0),g(E,m2,UU(1)),x1(1)-x1(2));
    Tline=@(E)reshape(T0(E),1,[]);

    En=linspace(0,Umax,NE);
    D=zeros(1,NE);
    for ii=1:NE
        T1=Tline(En(ii));
        D(ii)=abs(g(En,m1,-V*e))/abs(g(En,m1,0))*mm(1)/mm(5)*abs((T1(4)*T1(1)-T1(2)*T1(3))/T1(4))^2;
    end

    Emax=En(islocalmax(D));
    Dmax=D(islocalmax(D));
    if(Dmax(1)>Dmax(2))
        Ennew=linspace(max(Emax(1)-0.05*e,0),Emax(1)+0.05*e,2*NE);
    else
        Ennew=linspace(max(Emax(2)-0.05*e,0),Emax(2)+0.05*e,2*NE);
    end
    Dnew=zeros(1,2*NE);
    for ii=1:2*NE
        T1=Tline(Ennew(ii));
        Dnew(ii)=abs(g(Ennew,m1,-V*e))/abs(g(Ennew,m1,0))*mm(1)/mm(5)*abs((T1(4)*T1(1)-T1(2)*T1(3))/T1(4))^2;
    end

    S=@(Ez)log((1+exp((mu*e-Ez)/(kT)))./(1+exp((mu*e-Ez-V*e)/(kT))));
    S=2*pi*m1*kT/(hbar^2)*S(Ennew);
    J(i)=e/hbar*2/(2*pi)^3*trapz(-S.*Dnew,Ennew);

    if isnan(J(i)) && i>1
        J(i)=J(i-1);
        sprintf('J is nan');
    end

    if(i==1) 
        J(i)=0;
        D1=D;
    end

    f=figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.5]);
    subplot(1,3,1)
%     x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
%     y=[0 0 UU UU(end) UU(end)];
    x=x1*1e9 ; y=UU;
%     Emax=En(islocalmax(D1));
    plot(x,y/e,'-b')
    x=x+dx1*1e9/2;
    hold on
    text(7, 0.5, num2str(J(i)))
%     if (~isempty(Emax))
%         plot([x(1) x(end)]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
%     end
    ylim([-Umax/e-.1, Umax/e])
    xlim([x(1)-0.1 x(end)+0.1])
    title({'Энергетический профиль', 'структуры'})
    xlabel('x,нм')
    ylabel('E, эВ')
    grid on


    subplot(1,3,2)
    plot(D,En/e,'k','LineWidth',1);
    hold on
    fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
    plot(Dnew,Ennew/e,'r--','LineWidth',2);
    plot(S,Ennew/e);
    plot(S.*Dnew,Ennew/e);
    plot(Dmax,Emax/e,'o')
    ylim([-Umax/e-.1, Umax/e])
    xlim([0 1])
    ylabel('E,эВ')
    xlabel('D')
    title({'Коэффициент' ,'прозрачности'})
    grid on
    box on
    set(gca,'FontSize',11,'fontWeight','bold')

    subplot(1,3,3)
    plot(log(D),En/e,'k','LineWidth',1);
    hold on
    plot(log(Dnew),Ennew/e,'r--','LineWidth',2);
    plot(log(Dmax),Emax/e,'o')
    ylim([-Umax/e-.1, Umax/e])
    ylabel('E,эВ')
    xlabel('D')
    title({'Коэффициент' ,'прозрачности'})
    grid on
    box on
    set(gca,'FontSize',11,'fontWeight','bold')

    writeVideo(v,getframe(gcf));
    close(f)
end
close all
close(v)

%% 
UU(1:N1)=0;
UU(x1>=2e-9 & x1<4.5e-9)=U0;
UU(x1>=9.2e-9 & x1<11.7e-9)=U0;

figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.5]);
subplot(1,3,1)
x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
y=[0 0 UU UU(end) UU(end)];
Emax=En(islocalmax(D1));
plot(x,y,'-b')
x=x+dx1*1e9/2;
% hold on
% if (~isempty(Emax))
%     plot([x(1) x(end)]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
% end
ylim([-Umax/e-.1, Umax/e])
xlim([x(1)-0.1 x(end)+0.1])
title({'Энергетический профиль', 'структуры'})
xlabel('x,нм')
ylabel('E, эВ')
grid on


subplot(1,3,2)
hold on
fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
plot(D1,En/e,'k','LineWidth',1);
ylim([-Umax/e-.1, Umax/e])
xlim([0 1])
ylabel('E,эВ')
xlabel('D')
title({'Коэффициент' ,'прозрачности'})
grid on
box on
set(gca,'FontSize',11,'fontWeight','bold')


subplot(1,3,3)
hold on
plot(Volt,J,'k','LineWidth',1);
ylim([-2e8, 2e8])
ylabel('J, А/м^2')
xlabel('V, В')
title('ВАХ')
grid minor
box on
set(gca,'FontSize',11,'fontWeight','bold')


% %% graphic
% 
% f=figure('Units','normalized','OuterPosition',[0.1 0.1 0.45 0.45]);
% subplot(1,3,1)
% x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
% y=[0 0 UU/e -V -V];
% Emax=En(islocalmax(D)); 
% plot(x,y,'-b')
% x=x+dx1*1e9/2;
% hold on
% plot(x,y,'_','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor',[0.95 0.95 0.95],'MarkerIndices',1:4:N1+4)
% if (~isempty(Emax))
% plot([b,a+b]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
% end
% % text((x1(end)+x1(1))/2, -Umax*0.9, sprintf('V = %0.2f В', V),"FontSize",12,"FontWeight","bold",'HorizontalAlignment','center')
% ylim([-Umax/e-.1, Umax/e])
% xlim([x(1)-0.1 x(end)+0.1])
% title({'Энергетический профиль', 'структуры'})
% xlabel('x,нм')
% ylabel('E, эВ')
% grid on
% set(gca,'FontSize',12,'fontWeight','bold')
% 
% subplot(1,3,2)
% hold on
% fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
% plot(D,En/e,'k','LineWidth',1);
% plot(Dav,En/e,'--k','LineWidth',1);
% ylim([-Umax/e-.1, Umax/e])
% xlim([0 1])
% ylabel('E,эВ')
% xlabel('D')
% title({'Коэффициент' ,'прозрачности'})
% grid on
% box on
% legend('Факт. масса', "Ср. масса", 'Location','south')
% set(gca,'FontSize',12,'fontWeight','bold')
% 
% 
% subplot(1,3,3)
% hold on
% plot(Volt,J*1e6,'k','LineWidth',1);
% plot(Volt,Jav*1e6,'--k','LineWidth',1);
% ylim([ylim*[1;0]-ylim*[-1; 1]/3, ylim*[0;1]])
% % xlim([0 1])
% ylabel('J, мкА')
% xlabel('V, В')
% title('ВАХ')
% grid on
% box on
% legend('Факт. масса', "Ср. масса", 'Location','south')
% set(gca,'FontSize',12,'fontWeight','bold')
% 
% fprintf('done')