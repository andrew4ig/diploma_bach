clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;

% B=(2.5:0.1:10);
Xval=0.04:0.04:1;

eps=zeros(1,length(Xval));
Efact=zeros(1,length(Xval));
Eaver=zeros(1,length(Xval));
zlim=[-6,10];


v = VideoWriter('EpsilonXm.mp4','MPEG-4');
v.FrameRate = 8;
v.Quality = 100;
open(v);
% 
% 
for ind=1:length(Xval)
%     b=B(ind)*1e-9;
    X=Xval(ind);

%     X=1;
    %initial paramets
    a=3e-9;                 %size of chanel
    b=5e-9;                 %size of barrier
    U0=e*X*0.74;            %height of barrier,  eV
    m1=0.067*m0;            %eff. mass in GaAs, kg
    m2=(0.067+0.083*0.3)*m0;  %eff. mass in AlGaAs(x), kg

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
    
    Umax=1;
    N=2500;
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

    Efact(ind)=Emaxi(1)/e;

    if(isempty(Emaxi))
        En=U0;
    else
        En=Emaxi(1);
    end
    
    figure ('Units','normalized','OuterPosition', [0.05 0.05 0.45 0.45]);
    subplot(1,2,2)
    hold on
    plot(D,E/e,'k')
    if(~isempty(Emini))
        plot(Dmin, Emini/e, 'kv','HandleVisibility','off')
    end
    if (~isempty(Emaxi))
        plot(Dmax, Emaxi/e, 'kv','HandleVisibility','off')
    end
    
    subplot(1,2,1)
    S=zeros(1,100);  k1=gm1(En); k2=gm2(En);
    Psi1=@(A,B,x)A*exp(1i*k1*x)+B*exp(-1i*k1*x);
    Psi2=@(A,B,x)A*exp(k2*x)+B*exp(-k2*x);
    Nx=100;
    x=linspace(zlim(1),zlim(2),Nx)*1e-9;
    T=Koef(En);
    for i=1:Nx
        if(x(i)<0)
            S(i)=Psi1(A1,T(1),x(i));
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi2(T(2),T(3),x(i));
        elseif (x(i)>=b)
            S(i)=Psi1(T(4),0,x(i));
        end
    end
    plot(x*1e9,En/e+abs(S).^2/20,'k','LineWidth',3,'HandleVisibility','off')
    hold on
    plot(x*1e9,En/e+real(S)/20,'k--','HandleVisibility','off')
    plot(x*1e9,En/e+imag(S)/20,'k.','HandleVisibility','off')    
    xnew=[zlim(1)*1e-9,   0,  0,  b,  b, zlim(2)*1e-9];
    Unew=[0     0   U0  U0  0   0]/e;
    mass=[m m mb mb m m];
    plot(xnew*1e9, Unew,'HandleVisibility','off')
    xlim(zlim)
    ylim([0 Umax])
    yline(En/e,'Color',[0.7 0.7 0.7],'HandleVisibility','off');
    xlabel('x,????')
    ylabel('E, ????')
    text(-5, 0.4, sprintf('z = %0.2f ????% ', b*1e9),"FontSize",12,"FontWeight","bold")

    yyaxis right
    plot(xnew*1e9,mass/m0,'-');
    ylabel('m/m_0')
    ylim([0 0.5])

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
    

    D=zeros(1,N);
    for i=1:N
        T=Koef(E(i));
        D(i)=abs(T(4))^2;
    end
    Emaxi=E(islocalmax(D));
    Dmax=D(islocalmax(D));
    Emini=E(islocalmin(D));
    Dmin=D(islocalmin(D));

    Eaver(ind)=Emaxi(1)/e;
    
    if(isempty(Emaxi))
        En=U0;
    else
        En=Emaxi(1);
    end
        
    subplot(1,2,2)
    hold on
    plot(D,E/e,'r')
    hold on
    if(~isempty(Emini))
        plot(Dmin, Emini/e, 'r^','HandleVisibility','off')
    end
    if (~isempty(Emaxi))
        plot(Dmax, Emaxi/e, 'r^','HandleVisibility','off')
    end
    grid on
    ylabel('E,????')
    xlabel('D')
    title('?????????????????????? ????????????????????????????')
    ylim([0 Umax])
    xlim([0 1])
    legend('m_{????????????}','m_{????}', 'Location','northwest')
    box on
    
    subplot(1,2,1)
    yyaxis left
    S=zeros(1,100);  k1=gm1(En); k2=gm2(En);
    Psi1=@(A,B,x)A*exp(1i*k1*x)+B*exp(-1i*k1*x);
    Psi2=@(A,B,x)A*exp(k2*x)+B*exp(-k2*x);
%     x=linspace(zlim(1),zlim(2),N)*1e-9;
    T=Koef(En);
    for i=1:Nx
        if(x(i)<0)
            S(i)=Psi1(A1,T(1),x(i));
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi2(T(2),T(3),x(i));
        elseif (x(i)>=b)
            S(i)=Psi1(T(4),0,x(i));
        end
    end
    plot(x*1e9,En/e+abs(S).^2/20,'r','LineWidth',3,'HandleVisibility','off')
    hold on
    plot(x*1e9,En/e+real(S)/20,'r--','HandleVisibility','off')
    plot(x*1e9,En/e+imag(S)/20,'r.','HandleVisibility','off')    
    
    xnew=[zlim(1)*1e-9,   0,  0,  b,  b, zlim(2)*1e-9];
    mass=[m m mb mb m m];
    yline(En/e,'Color',[0.7 0.7 0.7],'HandleVisibility','off');

    yyaxis right
    plot(xnew*1e9,mass/m0,'--');

    legend('m_{????????????}/m_0->','m_{????}/m_0->')
    grid on
    title('???????????????? ???? ?? ??????????????')
    
    eps(ind)=abs((Efact(ind)-Eaver(ind))./Efact(ind)*100);

    subplot(1,2,2)
    text(0.4, 0.4, sprintf('eps = %0.1f%', eps(ind)),"FontSize",12,"FontWeight","bold")
    text(0.4, 0.3, sprintf('E_{fact} = %0.2f ????', Efact(ind)),"FontSize",12,"FontWeight","bold")
    yline(Efact(ind),'HandleVisibility','off')
    text(0.4, 0.2, sprintf('E_{aver} = %0.2f ????', Eaver(ind)),"FontSize",12,"FontWeight","bold")
    yline(Eaver(ind),'HandleVisibility','off')


    s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
    s.matlab.fonts.editor.normal.Size.TemporaryValue = 12;

    writeVideo(v,getframe(gcf));

    close all
end
% close all
close(v)

fprintf('done')

%% graph
figure ('Units','normalized','OuterPosition', [0.05 0.05 0.45 0.45]);
plot(Xval,eps,'k')
grid on
ylabel('eps, %')
xlabel('x')
title({'?????????????????????? ???????????????????? ?? ??????????????????????', '???? ???????????? ?????????????? ?????? b=5 ????','?????? ?????????????????????????? ???????? ?????????????????????? ?????????????????? x=0.3'})

yyaxis right
plot(Xval,Efact,'-');
hold on
plot(Xval,Eaver,'--');
plot(Xval,abs(Eaver-Efact),'-.')
legend('??????????????.','????????????. ??????.', '??????????. ??????.', 'dE_{res}')
ylabel('E_{res}, ????')
box on
xlim([Xval(1), Xval(end)])
  

s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
s.matlab.fonts.editor.normal.Size.TemporaryValue = 14;
