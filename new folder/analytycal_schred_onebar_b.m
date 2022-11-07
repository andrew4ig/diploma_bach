clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;

B=3:0.1:10;

v = VideoWriter('OneBarB.mp4','MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v);

Ex=zeros(1,length(B));
Wx=zeros(1,length(B));
Dx=zeros(1,length(B));

for ind=1:length(B)
    b=B(ind)*1e-9;

    X=0.3;
    %initial paramets
    a=3e-9;                 %size of chanel, nm ?????????????
    % b=2e-9;                 %size of barrier, nm ?????????????
    U0=e*X*0.74;            %height of barrier,  eV
    m1=0.067*m0;            %eff. mass in GaAs, kg
    m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg

    %solving Schr eq
    m=m1; mb=m1;
    m1=(m1+m2)/2; m2=m1;
    gm1=@(E)sqrt(2*m*(E))/hbar;
    gm2=@(E)sqrt(2*mb*(-E+U0))/hbar;
    gm3=@(E)sqrt(2*m*(E))/hbar;
    % gm4=@(E)sqrt(2*mb*(-E+U0))/hbar;
    % gm5=@(E)sqrt(2*m*(E))/hbar;
    A1=1;

    MKoef=@(E)[ 1,              -1,                         -1,                                 0;
        -1i*gm1(E)/m,   -gm2(E)/mb,                 gm2(E)/mb,                          0;
        0,              exp(gm2(E)*b),              exp(-gm2(E)*b),                     -exp(1i*gm3(E)*b);
        0,              gm2(E)/mb.*exp(gm2(E)*b),   -gm2(E)/mb.*exp(-gm2(E)*b),         -1i*gm3(E)/m.*exp(1i*gm3(E)*b);];
    MSvob=@(E)[-A1; -A1*1i*gm1(E)/m; 0; 0];
    % Koef=@(E)bicgstab(MKoef(E),MSvob(E), 1e-12, 100);
    % Koef=@(E)gmres(MKoef(E),MSvob(E),10,1e-12,15);
    Koef=@(E)(MKoef(E)\MSvob(E));

    Umax=2;
    N=300;
    E=linspace(0,Umax*e,N);
    D=zeros(1,N);
    for i=1:N
        T=Koef(E(i));
        D(i)=abs(T(4))^2;
    end
    E0=E(islocalmax(D));
    D0=D(islocalmax(D));

    Emax=hbar^2/(2*mb)*(pi*(1:3)/b).^2+U0;
    Emin=hbar^2/(2*mb)*(pi/b).^2.*((1:3)+1/2).^2+U0;
    DbyE=@(E) abs((4*gm1(E).*1i.*gm2(E).*exp(-1i*gm1(E)*b))./((gm1(E)+1i*gm2(E)).^2.*exp(-1i*1i*gm2(E)*b)-(gm1(E)-1i*gm2(E)).^2.*exp(1i*1i*gm2(E)*b))).^2;
    

    Emax=Emax(Emax<Umax*e);
    Emin=Emin(Emin<Umax*e);

    Eupper=linspace(Emax(1)*0.5,Emin(1)+0.01*e,N);
    Dupper=zeros(1,N);
    for i=1:N
        T=Koef(Eupper(i));
        Dupper(i)=abs(T(4))^2;
    end

    if(isempty(E0))
        En=U0/2;
    else
        En=E0(1);
    end

    f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
    subplot(1,2,2)
    hold on
    plot(D,E/e,'k')
    hold on
    %     plot(D(id1),E(id1)/e,'--r','LineWidth',2)
    %     plot(D(id2),E(id2)/e,'--b','LineWidth',2)
    width=E(islocalmin(abs(D-1+exp(-1))));
    %     width=[width(1),width(2)]/e;
    %     plot((1-exp(-1))*[1,1],width,'k','LineWidth',2)
    %     label = {sprintf(' E=%.0fмэВ \n D=%.1f \n g=%.1fмэВ',[E0(1)/e*1000,D0(1),(width(2)-width(1))*1000]),sprintf('E = %.0f мэВ\nD=%.1f',[E0(2)/e*1000,D0(2)])};
    %     yl=yline(E0/e,'--',label,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top');
%     if(~isempty(Emin))
%         yline(Emin/e,'r--', 'Min', 'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
%         plot(DbyE(Emin), Emin/e, 'r^')
%     end
%     if (~isempty(Emax))
%         yline(Emax/e,'b--', 'Max', 'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left')
%         plot(DbyE(Emax), Emax/e, 'bv')
%     end
    grid on
%     plot(DbyE(E), E/e, '--g')
    ylabel('E,эВ')
    xlabel('D')
    title('Коэффициент прохождения')
    ylim([0 1.5])
    xlim([0 1])
    box on

    subplot(1,2,1)
    S=zeros(1,100);  k1=gm1(En); k2=gm2(En);
    Psi1=@(A,B,x)A*exp(1i*k1*x)+B*exp(-1i*k1*x);
    Psi2=@(A,B,x)A*exp(k2*x)+B*exp(-k2*x);
    N=100;
    x=linspace(-5e-9,b+5e-9,N);
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
    xnew=[x(1),   0,  0,  b,  b, x(end)];
    Unew=[0     0   U0  U0  0   0]/e;
    plot(xnew*1e9, Unew,'HandleVisibility','off')
    ylim([0 1.5])
    xlim([x(1) x(end)]*1e9)
    Eline=yline(En/e,'Color',[0.7 0.7 0.7]);
    xlabel('x,нм')
    ylabel('E, эВ')
    legend('|\Psi|^2','Re(\Psi)','Im(\Psi)','E_0')
    grid on
    title('Динамика ВФ в профиле')

    if ind==0
        Wx(ind)=0;
    else
        Dmintemp=Dupper(islocalmin(Dupper));
        Dmaxtemp=Dupper(islocalmax(Dupper));
        Dtemp=Dmintemp(1)+(Dmaxtemp(1)-Dmintemp(1))*(1-exp(-1));
        Dtemp=abs(Dupper-Dtemp);
        Etemp=Eupper(islocalmin(Dtemp));
       

        subplot(1,2,2)
        hold on
        xline((Dmintemp(1)+(Dmaxtemp(1)-Dmintemp(1))*(1-exp(-1))))
        plot((Dmintemp(1)+(Dmaxtemp(1)-Dmintemp(1))*(1-exp(-1)))*[1 1],[Etemp(2) Etemp(1)]/e,'o')

        Etemp=Etemp(2)-Etemp(1);
        Wx(ind)=Etemp;
    end

    Ex(ind)=Eupper(1);
    Dupper=Dupper(islocalmax(Dupper));
    Dx(ind)=Dupper(1);

    writeVideo(v,getframe(gcf));
    close(f)
end
close all
close(v)

%% peaks widths ans values
f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
plot(B, Ex/e,'k')
hold on
plot(B,Wx/e,'k--')
ylabel('E, эВ')

yyaxis right
plot(B, Dx)
ylabel('D')
ylim([0 1.5])

xlim([B(1), B(end)])
box on
grid on
legend('E_{рез}','Уширение','К-т пр-ти')
xlabel('b, нм')
title({'Завсимость резонасного коэффициента прохождения', 'от ширины барьера'})
s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
s.matlab.fonts.editor.normal.Size.TemporaryValue = 14;
