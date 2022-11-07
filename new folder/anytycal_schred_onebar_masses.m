clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;


% v = VideoWriter('OneBarMasses.mp4','MPEG-4');
% v.FrameRate = 8;
% v.Quality = 100;
% open(v);

zlim=[-6,10];

% for b=[(4:0.05:7)*1e-9]

    X=0.3;
    %initial paramets
    a=3e-9;                 %size of chanel, nm ?????????????
    b=8e-9;                 %size of barrier, nm ?????????????
    U0=e*1;%*X*0.74;            %height of barrier,  eV
    m1=0.067*m0;            %eff. mass in GaAs, kg
    m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg
    
    %solving Schr eq
    m=m1; mb=m2;
    gm1=@(E)sqrt(2*m*(E))/hbar;
    gm2=@(E)sqrt(2*mb*(-E+U0))/hbar;
    gm3=@(E)sqrt(2*m*(E))/hbar;
    A1=1;
    
    MKoef=@(E)[ 1,              -1,                         -1,                                 0;
                -1i*gm1(E)/m,   -gm2(E)/mb,                 gm2(E)/mb,                          0;
                0,              exp(gm2(E)*b),              exp(-gm2(E)*b),                     -exp(1i*gm3(E)*b);
                0,              gm2(E)/mb.*exp(gm2(E)*b),   -gm2(E)/mb.*exp(-gm2(E)*b),         -1i*gm3(E)/m.*exp(1i*gm3(E)*b);];
    MSvob=@(E)[-A1; -A1*1i*gm1(E)/m; 0; 0];
    Koef=@(E)(MKoef(E)\MSvob(E));
    
    Umax=3;
    N=300;
    E=linspace(0,Umax*e,N);
    D=zeros(1,N);
    for i=1:N
        T=Koef(E(i));
        D(i)=abs(T(4))^2;
    end
    
    Emaxi=E(islocalmax(D));
    Dmax=D(islocalmax(D));
    Emini=E(islocalmin(D));
    Dmin=D(islocalmin(D));

    Efact=Emaxi(1);

    if(isempty(Emaxi))
        En=U0;
    else
        En=Emaxi(1);
    end
    
    figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.9]);
    subplot(2,2,2)
    hold on
    plot(D,E/e,'k')
    hold on
    if(~isempty(Emini))
        yline(Emini/e,'r--', 'Min', 'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
        plot(Dmin, Emini/e, 'r^')
    end
    if (~isempty(Emaxi))
        yline(Emaxi/e,'b--', 'Max', 'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
        plot(Dmax, Emaxi/e, 'bv')
    end
    grid on
    ylabel('E,эВ')
    xlabel('D')
    title('Коэффициент прохождения')
    ylim([0 3])
    xlim([0 1])
    box on
    
    subplot(2,2,1)
    S=zeros(1,100);  k1=gm1(En); k2=gm2(En);
    Psi1=@(A,B,x)A*exp(1i*k1*x)+B*exp(-1i*k1*x);
    Psi2=@(A,B,x)A*exp(k2*x)+B*exp(-k2*x);
    N=100;
    x=linspace(zlim(1),zlim(2),N)*1e-9;
    T=Koef(En);
    for i=1:N
        if(x(i)<0)
            S(i)=Psi1(A1,T(1),x(i));
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi2(T(2),T(3),x(i));
        elseif (x(i)>=b)
            S(i)=Psi1(T(4),0,x(i));
        end
    end
    plot(x*1e9,En/e+abs(S).^2/40,'k','LineWidth',3)
    hold on
    plot(x*1e9,En/e+real(S)/40,'b--')
    plot(x*1e9,En/e+imag(S)/40,'r-.')    
    xnew=[zlim(1)*1e-9,   0,  0,  b,  b, zlim(2)*1e-9];
    Unew=[0     0   U0  U0  0   0]/e;
    mass=[m m mb mb m m];
    plot(xnew*1e9, Unew,'HandleVisibility','off')
    xlim(zlim)
    ylim([0 3])
    yline(En/e,'Color',[0.7 0.7 0.7]);
    xlabel('x,нм')
    ylabel('E, эВ')
    text(-5, 0.7, sprintf('b = %0.2f нм', b*1e9),"FontSize",12,"FontWeight","bold")

    yyaxis right
    plot(xnew*1e9,mass/m0);
    ylabel('m/m_0')
    ylim([0 0.5])

    legend('|\Psi|^2','Re(\Psi)','Im(\Psi)','E_0', 'm/m0')
    grid on
    title('Динамика ВФ в профиле')
    

    %solving Schr eq
    m=(m1+m2)/2; mb=m;
    gm1=@(E)sqrt(2*m*(E))/hbar;
    gm2=@(E)sqrt(2*mb*(-E+U0))/hbar;
    gm3=@(E)sqrt(2*m*(E))/hbar;
    A1=1;
    
    MKoef=@(E)[ 1,              -1,                         -1,                                 0;
                -1i*gm1(E)/m,   -gm2(E)/mb,                 gm2(E)/mb,                          0;
                0,              exp(gm2(E)*b),              exp(-gm2(E)*b),                     -exp(1i*gm3(E)*b);
                0,              gm2(E)/mb.*exp(gm2(E)*b),   -gm2(E)/mb.*exp(-gm2(E)*b),         -1i*gm3(E)/m.*exp(1i*gm3(E)*b);];
    MSvob=@(E)[-A1; -A1*1i*gm1(E)/m; 0; 0];
    Koef=@(E)(MKoef(E)\MSvob(E));
    
    Umax=3;
    N=300;
    E=linspace(0,Umax*e,N);
    D=zeros(1,N);
    for i=1:N
        T=Koef(E(i));
        D(i)=abs(T(4))^2;
    end
    Emaxi=E(islocalmax(D));
    Dmax=D(islocalmax(D));
    Emini=E(islocalmin(D));
    Dmin=D(islocalmin(D));

    Eav=Emaxi(1);
    
    if(isempty(Emaxi))
        En=U0;
    else
        En=Emaxi(1);
    end
        
    subplot(2,2,4)
    hold on
    plot(D,E/e,'k')
    hold on
    if(~isempty(Emini))
        yline(Emini/e,'r--', 'Min', 'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
        plot(Dmin, Emini/e, 'r^')
    end
    if (~isempty(Emaxi))
        yline(Emaxi/e,'b--', 'Max', 'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
        plot(Dmax, Emaxi/e, 'bv')
    end
    grid on
    ylabel('E,эВ')
    xlabel('D')
    title('Коэффициент прохождения')
    ylim([0 3])
    xlim([0 1])
    box on
    
    subplot(2,2,3)
    S=zeros(1,100);  k1=gm1(En); k2=gm2(En);
    Psi1=@(A,B,x)A*exp(1i*k1*x)+B*exp(-1i*k1*x);
    Psi2=@(A,B,x)A*exp(k2*x)+B*exp(-k2*x);
    N=100;
    x=linspace(zlim(1),zlim(2),N)*1e-9;
    T=Koef(En);
    for i=1:N
        if(x(i)<0)
            S(i)=Psi1(A1,T(1),x(i));
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi2(T(2),T(3),x(i));
        elseif (x(i)>=b)
            S(i)=Psi1(T(4),0,x(i));
        end
    end
    plot(x*1e9,En/e+abs(S).^2/40,'k','LineWidth',3)
    hold on
    plot(x*1e9,En/e+real(S)/40,'b--')
    plot(x*1e9,En/e+imag(S)/40,'r-.')    
    
    xnew=[zlim(1)*1e-9,   0,  0,  b,  b, zlim(2)*1e-9];
    Unew=[0     0   U0  U0  0   0]/e;
    mass=[m m mb mb m m];
    plot(xnew*1e9, Unew,'HandleVisibility','off')
    xlim(zlim)
    ylim([0 3])
    yline(En/e,'Color',[0.7 0.7 0.7]);
    xlabel('x,нм')
    ylabel('E, эВ')
    text(-5, 0.7, sprintf('b = %0.2f нм', b*1e9),"FontSize",12,"FontWeight","bold")

    yyaxis right
    plot(xnew*1e9,mass/m0);
    ylabel('m/m_0')
    ylim([0 0.5])

    legend('|\Psi|^2','Re(\Psi)','Im(\Psi)','E_0', 'm/m0')
    grid on
    title('Динамика ВФ в профиле')
    
%     
%     writeVideo(v,getframe(gcf));
% end

% close all
% close(v)
fprintf('done')

Efact/e
Eav/e
(Efact-Eav)/Efact*100
