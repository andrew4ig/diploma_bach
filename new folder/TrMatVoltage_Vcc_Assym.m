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

% X1=1;
% X2=1;
a=3e-9;
b1=2e-9;
b2=2e-9;
% U1=X1*0.74;
% U2=X2*0.74;
m1=0.067*m0;
% m21=(0.067+0.083*X1)*m0;
% m22=(0.067+0.083*X2)*m0;

dx1=2e-10;
L1=b1+a+b2;
NE=1500;
N1=floor(L1/dx1);
x1=(0:N1-1)*dx1;

Vmax=0.8;
Volt=0:0.02:Vmax;
Umax=1*e;
En=linspace(0, Umax,NE); 
dE=(En(2)-En(1));
n=1e21;
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.3;



% [J, D]=current(mm, m1, m21, m22, Volt, En, dE, NE, a, b1, b2, U1, U2, e, kT, hbar, mu, N1, x1);

% v = VideoWriter('MassesAproach.mp4','MPEG-4');
% v.FrameRate = 4;
% v.Quality = 100;
% open(v); 
Xval=1;i=1;

    %% calculate 1
    X1=Xval(i);
    X2=0.3;
    U1=X1*0.74;
    U2=X2*0.74;
    m21=(0.067+0.083*X1)*m0;
    m22=(0.067+0.083*X2)*m0;
    mm=ones(1,N1)*m1;
    mm(x1<b1)=m21;
    mm(x1>b1+a)=m22;
    [Jleft, Dleft]=current(mm, m1, m21, m22, Volt, En, dE, NE, a, b1, b2, U1, U2, e, kT, hbar, mu, N1, x1);
    
    
    %% graph1
    UU=zeros(1,N1);
    UU(x1<b1)=U1;
    UU(x1>b1+a)=U2;
    
    f=figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.5]);
    subplot(2,6,1:2)
    x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
    y=[0 0 UU 0 -0];
    Emax=En(islocalmax(Dleft));
    plot(x,y,'-b')
    x=x+dx1*1e9/2;
    hold on
    if (~isempty(Emax))
        plot([b1,a+b1]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
    end
    ylim([-Umax/e-.1, Umax/e])
    xlim([x(1)-0.1 x(end)+0.1])
    title({'Энергетический профиль', 'структуры'})
    xlabel('x,нм')
    ylabel('E, эВ')
    grid on
    
    subplot(2,6,3:4)
    hold on
    % fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
    plot(Dleft,En/e,'b','LineWidth',1);
    ylim([-Umax/e-.1, Umax/e])
    xlim([0 1])
    ylabel('E,эВ')
    xlabel('D')
    title({'Коэффициент' ,'прозрачности'})
    grid on
    box on
    
    subplot(2,6,7:9)
    hold on
    plot(Volt,Jleft,'b','LineWidth',1);
    ylabel('J, А/м^2')
    xlabel('V, В')
    title('ВАХ')
    ylim([0, 1e11])
    grid on
    box on

    %% calculate 2
    X1=0.3;
    X2=Xval(i);
    U1=X1*0.74;
    U2=X2*0.74;
    m21=(0.067+0.083*X1)*m0;
    m22=(0.067+0.083*X2)*m0;
    mm=ones(1,N1)*m1;
    mm(x1<b1)=m21;
    mm(x1>b1+a)=m22;
    mm=ones(1,N1)*m1;
    mm(x1<b1)=m21;
    mm(x1>b1+a)=m22;
    [Jright, Dright]=current(mm, m1, m21, m22, Volt, En, dE, NE, a, b1, b2, U1, U2, e, kT, hbar, mu, N1, x1);
    
    %% graph2
    UU=zeros(1,N1);
    UU(x1<b1)=U1;
    UU(x1>b1+a)=U2;
    
    subplot(2,6,5:6)
    x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
    y=[0 0 UU 0 -0];
    Emax=En(islocalmax(Dright));
    plot(x,y,'.-r','LineWidth',2)
    x=x+dx1*1e9/2;
    hold on
    if (~isempty(Emax))
        plot([b1,a+b1]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
    end
    ylim([-Umax/e-.1, Umax/e])
    xlim([x(1)-0.1 x(end)+0.1])
    title({'Энергетический профиль', 'структуры'})
    xlabel('x,нм')
    ylabel('E, эВ')
    grid on
    
    subplot(2,6,3:4)
    hold on
    % fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
    plot(Dright,En/e,'-.r','LineWidth',2);
    
    subplot(2,6,10:12)
    hold on
    plot(Volt,Jright,'r','LineWidth',1);
    ylabel('J, А/м^2')
    xlabel('V, В')
    title('ВАХ')
    ylim([0, 1e11])
    grid on
    box on

%     writeVideo(v,getframe(gcf));
%     close(f)
% end
% close(v)



%% function calculating transparency and current
function [J,D1] = current(mm,m1,m21,m22, Volt, En,dE, NE, a, b1, b2, U1, U2, e, kT, hbar, mu, N1, x1)
    g=@(E,m,U)sqrt(2*m*(E-U))/hbar;
    T=@(m1,m2,g1,g2,z)[0.5*(1+g1/g2*m2/m1)*exp(-1i*(g2-g1)*z),   0.5*(1-g1/g2*m2/m1)*exp(-1i*(g2+g1)*z);
                       0.5*(1-g1/g2*m2/m1)*exp( 1i*(g2+g1)*z),   0.5*(1+g1/g2*m2/m1)*exp( 1i*(g2-g1)*z)];
    J=zeros(1,length(Volt));
%     D1=zeros(1,NE);

    for i=1:length(Volt)
        V=Volt(i);

        UU=zeros(1,N1);
        
        for ii=1:N1
            if(x1(ii)>=0&&x1(ii)<b1)
                UU(ii)=U1-V*ii/N1;
            elseif (x1(ii)>=b1&&x1(ii)<b1+a)
                UU(ii)=-V*ii/N1;
            elseif(x1(ii)>=b1+a)
                UU(ii)=U2-V*ii/N1;
            end
        end
        UU=UU*e;
%         UU=(UU-V/N1/2)*e;
        
        T0=@(E)T(m22,m1,g(E,m22,UU(N1)),g(E,m1,-V*e),x1(N1));
        for ii=N1-1:-1:1
            T0=@(E)T0(E)*T(mm(ii),mm(ii+1),g(E,mm(ii),UU(ii)),g(E,mm(ii+1),UU(ii+1)),x1(ii));
        end
        T0=@(E)T0(E)*T(m1,m21,g(E,m1,0),g(E,m21,UU(1)),x1(1)-x1(2));
        Tline=@(E)reshape(T0(E),1,[]);
    
        D=zeros(1,NE);
        for ii=1:NE
            T1=Tline(En(ii));
            D(ii)=abs(g(En,m1,-V*e))/abs(g(En,m1,0))*mm(1)/mm(5)*abs((T1(4)*T1(1)-T1(2)*T1(3))/T1(4))^2;
        end
       
        S=@(Ez)log((1+exp((mu*e-Ez)/(kT)))./(1+exp((mu*e-Ez-V*e)/(kT))));
        S=2*pi*m1*kT/(hbar^2)*S(En);
        J(i)=e/hbar*2/(2*pi)^3*trapz(S.*D)*dE; 
    
        if isnan(J(i))
            if(i>1)
                J(i)=J(i-1);
            elseif(i==1)
                J(i)=0;
            end
        end

        if(i==1)
            D1=D;
        end
    end
end

