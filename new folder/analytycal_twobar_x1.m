clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;

Xval=0:0.02:1;

v = VideoWriter('TwoBarX1.mp4','MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v);

Ex=zeros(1,length(Xval));
Wx=zeros(1,length(Xval));
Dx=zeros(1,length(Xval));

for ind=1:length(Xval)
    X=Xval(ind);
    X2=0.3;

    %initial paramets
    a=3e-9;                 %size of chanel, nm ?????????????
    b=2e-9;                 %size of barrier, nm ?????????????
    b2=b;
    U0=e*X*0.74;            %height of barrier,  eV
    m1=0.067*m0;            %eff. mass in GaAs, kg
    m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg
    
    
    %solving Schr eq
    m=m1; mb=m2;
    U1=U0; 
    U2=e*X2*0.74;
    m2b=(0.067+0.083*X2)*m0;
    gm1=@(E)sqrt(2*m*(E))/hbar;
    gm2=@(E)sqrt(2*mb*(-E+U1))/hbar;
    gm3=@(E)sqrt(2*m*(E))/hbar;
    gm4=@(E)sqrt(2*m2b*(-E+U2))/hbar;
    gm5=@(E)sqrt(2*m*(E))/hbar;
    A1=1;
    
    MKoef=@(E)[ 1,              -1,                         -1,                                 0,                                  0,                                  0,                                   0,                                   0;
                -1i*gm1(E)/m,   -gm2(E)/mb,                 gm2(E)/mb,                          0,                                  0,                                  0,                                   0,                                   0;
                0,              exp(gm2(E)*b),              exp(-gm2(E)*b),                     -exp(1i*gm3(E)*b),                  -exp(-1i*gm3(E)*b),                 0,                                   0,                                   0;
                0,              gm2(E)/mb.*exp(gm2(E)*b),   -gm2(E)/mb.*exp(-gm2(E)*b),         -1i*gm3(E)/m.*exp(1i*gm3(E)*b),     1i*gm3(E)/m.*exp(-1i*gm3(E)*b),     0,                                   0,                                   0;
                0,              0,                          0,                                  exp(1i*gm3(E)*(a+b)),               exp(-1i*gm3(E)*(a+b)),              -exp(gm4(E)*(a+b)),                  -exp(-gm4(E)*(a+b)),                 0;
                0,              0,                          0,                                  1i*gm3(E)/m.*exp(1i*gm3(E)*(a+b)), -1i*gm3(E)/m.*exp(-1i*gm3(E)*(a+b)), -gm4(E)/m2b.*exp(gm4(E)*(a+b)),       gm4(E)/m2b.*exp(-gm4(E)*(a+b)),       0;
                0,              0,                          0,                                  0,                                  0,                                  exp(gm4(E)*(a+b+b2)),                 exp(-gm4(E)*(a+b+b2)),                -exp(1i*gm5(E)*(a+b+b2));
                0,              0,                          0,                                  0,                                  0,                                  gm4(E)/m2b.*exp(gm4(E)*(a+b+b2)),      -gm4(E)/m2b.*exp(-gm4(E)*(a+b+b2)),    -1i*gm5(E)/m.*exp(1i*gm5(E)*(a+b+b2))];
    MSvob=@(E)[-A1; -A1*1i*gm1(E)/m; 0; 0; 0; 0; 0; 0];
    % Koef=@(E)bicgstab(MKoef(E),MSvob(E), 1e-12, 100);
    % Koef=@(E)gmres(MKoef(E),MSvob(E),10,1e-12,15);
    Koef=@(E)(MKoef(E)\MSvob(E));
    
    Umax=2.6;%max(U1,U2)*2/e;
    N=2000;
    E=linspace(0,Umax*e,N);
    D=zeros(1,N);
    for i=1:N
        T=Koef(E(i));
        D(i)=abs(T(8))^2;
    end
    Emax=E(islocalmax(D));
    Dmax=D(islocalmax(D));
    Emin=E(islocalmin(D));
    
%     En=Emax(1);
    Eupper=linspace(Emax(1)*0.5,Emin(1)+0.01*e,N);
    Dupper=zeros(1,N);
    for i=1:N
        T=Koef(Eupper(i));
        Dupper(i)=abs(T(8))^2;
    end
  
    f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
    subplot(1,2,2)
    plot(D,E/e,'b')
    hold on
    plot(Dupper,Eupper/e,'b')
%     plot(D(id1),E(id1)/e,'--r','LineWidth',2)
%     plot(D(id2),E(id2)/e,'--b','LineWidth',2)
%     width=E(islocalmin(abs(D-1+exp(-1))));
%     width=[width(1),width(2)]/e;
%     plot((1-exp(-1))*[1,1],width,'k','LineWidth',2)
%     label = {sprintf(' E=%.0fмэВ \n D=%.1f \n g=%.1fмэВ',[E0(1)/e*1000,D0(1),(width(2)-width(1))*1000]),sprintf('E = %.0f мэВ\nD=%.1f',[E0(2)/e*1000,D0(2)])};
%     yline(Efact(1)/e,'--',label,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top');
    grid on
    ylabel('E,эВ')
    xlabel('D')
    title('Коэффициент прохождения')
    ylim([0 Umax])
    xlim([0 1])
    hold off

%     axes('Position',[0.75 0.35 0.15 0.2])
%     box on
%     indexs = D>0.5 & E/e<0.3;
%     plot(D(indexs),E(indexs)/e,'k')
%     hold on
%     plot(D(indexs & E<E(islocalmin(D))),E(indexs & E<E(islocalmin(D)))/e,'--r','LineWidth',2)
%     plot((1-exp(-1))*[1,1],width,'k','LineWidth',2)
%     grid on
%     axis tight

    
    En=Eupper(islocalmax(Dupper));
    En=En(1);    
    subplot(1,2,1)
    S=zeros(1,100);
    k1=gm1(En); k2=gm2(En); k3=gm3(En); k4=gm4(En); k5=gm5(En);
    Psi=@(A,B,x,k)A*exp(k*x)+B*exp(-k*x);
    Psi2=@(A,B,x,k)A*exp(1i*k*x)+B*exp(-1i*k*x);
    Nx=300;
    %x=linspace(-b,a+3*b,Nx);
    x=linspace(-5e-9,b+b+a+5e-9,Nx);
%     x=linspace(-5e-9,25e-9,Nx);%b+b+a+5e-9
    T=Koef(En);
    % T=T./sqrt(T.*conj(T));
    for i=1:Nx
        if(x(i)<0)
            S(i)=Psi2(A1,T(1),x(i),k1);
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi(T(2),T(3),x(i),k2);
        elseif (x(i)>=b && x(i)<b+a)
            S(i)=Psi2(T(4),T(5),x(i),k3);
        elseif (x(i)>=b+a && x(i)<a+b+b2)
            S(i)=Psi(T(6),T(7),x(i),k4);
        elseif (x(i)>=a+b+b2)
            S(i)=Psi2(T(8),0,x(i),k5);
        end
    end
    plot(x*1e9,En/e+abs(S).^2/100,'k','LineWidth',2)
    hold on
    plot(x*1e9,En/e+real(S)/100,'b--')
    plot(x*1e9,En/e+imag(S)/100,'r-.')    
    xnew=[-b,   0,  0,  b,  b, a+b, a+b,    a+b+b2, a+b+b2, a+2*b+b2];
    Unew=[0     0   U1  U1  0   0   U2      U2      0       0]/e;
    plot(xnew*1e9, Unew,'HandleVisibility','off')
    xlim([x(1) x(end)]*1e9)
    ylim([0 Umax])
    yline(En/e,'Color',[0.7 0.7 0.7]);
    xlabel('x,нм')
    ylabel('E, эВ')
    legend('|\Psi|^2','Re(\Psi)','Im(\Psi)','E_0')
    grid on
    title('Динамика ВФ в профиле')
    
    if ind==1
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
    close (f)
end
close (v)

Xval(1)=[];
Dx(1)=[];
Ex(1)=[];
Wx(1)=[];

%% peaks widths ans values


f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
plot(Xval, Ex/e,'k')
hold on
plot(Xval,Wx/e,'k--')
ylabel('E, эВ')

yyaxis right
plot(Xval, Dx)
ylabel('E_{res}, эВ')
ylim([0 1.5])

xlim([Xval(1), Xval(end)])
box on
grid on
legend('E_{рез}','Уширение','К-т пр-ти')
xlabel('X')
title({'Зависмость резонасного коэффициента прохождения', 'от доли замещающего материала первого барьера'})
s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
s.matlab.fonts.editor.normal.Size.TemporaryValue = 14;