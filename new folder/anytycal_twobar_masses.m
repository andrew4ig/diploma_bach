clc
clear
close all



A=2:0.1:10;
Efct=zeros(1,length(A));
Eaver=zeros(1,length(A));

v = VideoWriter('TwoBarMassesA.mp4','MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v);



for ind=1:length(A)
    a=A(ind)*1e-9;
    %constants
    k=1.38e-23;
    hbar=1.0546e-34;
    m0=9.1e-31;
    e=1.6e-19;
    eps0=8.85e-12;

    %initial parametrs
    X=0.3;
    a=3e-9;                 %size of chanel, nm
    b=5e-9;                 %size of barrier, nm
    U0=e*X*0.74;            %height of barrier,  eV
    m1=0.067*m0;            %eff. mass in GaAs, kg
    m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg
    
    %solving Schr eq
    m=m1; mb=m2;
    U1=U0; 
    U2=U0;
    gm1=@(E)sqrt(2*m*(E))/hbar;
    gm2=@(E)sqrt(2*mb*(-E+U1))/hbar;
    gm3=@(E)sqrt(2*m*(E))/hbar;
    gm4=@(E)sqrt(2*mb*(-E+U2))/hbar;
    gm5=@(E)sqrt(2*m*(E))/hbar;
    A1=1;
    
    Z=@(E)[     1,              -1,                         -1,                                 0,...
                0,                                  0,                                   0,                                   0;
                -1i*gm1(E)/m,   -gm2(E)/mb,                 gm2(E)/mb,                          0,...
                0,                                  0,                                   0,                                   0;
                0,              exp(gm2(E)*b),              exp(-gm2(E)*b),                     -exp(1i*gm3(E)*b),...
                -exp(-1i*gm3(E)*b),                 0,                                   0,                                   0;
                0,              gm2(E)/mb.*exp(gm2(E)*b),   -gm2(E)/mb.*exp(-gm2(E)*b),         -1i*gm3(E)/m.*exp(1i*gm3(E)*b),     1i*gm3(E)/m.*exp(-1i*gm3(E)*b),...
                0,                                   0,                                   0;
                0,              0,                          0,                                  exp(1i*gm3(E)*(a+b)),               exp(-1i*gm3(E)*(a+b)),...
                -exp(gm4(E)*(a+b)),                  -exp(-gm4(E)*(a+b)),                 0;
                0,              0,                          0,                                  1i*gm3(E)/m.*exp(1i*gm3(E)*(a+b)), -1i*gm3(E)/m.*exp(-1i*gm3(E)*(a+b)),...
                -gm4(E)/mb.*exp(gm4(E)*(a+b)),       gm4(E)/mb.*exp(-gm4(E)*(a+b)),       0;
                0,              0,                          0,                                  0,                                  0,...
                exp(gm4(E)*(a+2*b)),                 exp(-gm4(E)*(a+2*b)),                -exp(1i*gm5(E)*(a+2*b));
                0,              0,                          0,                                  0,                                  0,...
                gm4(E)/mb.*exp(gm4(E)*(a+2*b)),      -gm4(E)/mb.*exp(-gm4(E)*(a+2*b)),    -1i*gm5(E)/m.*exp(1i*gm5(E)*(a+2*b))];
    Y=@(E)[-A1; -A1*1i*gm1(E)/m; 0; 0; 0; 0; 0; 0];
    X=@(E)(Z(E)\Y(E));
    
    N=3000;
    E=linspace(0,U0*e,N);
    D=zeros(1,N);
    for i=1:N
        T=X(E(i));
        D(i)=abs(T(8))^2;
    end


    Efct(ind)=Efact(1);

      
    f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);
    subplot(1,2,2)
    hold on
    plot(Dfact,E/e,'b')
    hold on
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

    En=Efact(1);
    subplot(1,2,1)
    S=zeros(1,100);
    k1=gm1(En); k2=gm2(En); k3=gm3(En); k4=gm4(En); k5=gm5(En);
    Psi=@(A,B,x,k)A*exp(k*x)+B*exp(-k*x);
    Psi2=@(A,B,x,k)A*exp(1i*k*x)+B*exp(-1i*k*x);
    Nx=300;
    %x=linspace(-b,a+3*b,Nx);
    x=linspace(-b,20e-9,Nx);
    T=Koef(En);
    % T=T./sqrt(T.*conj(T));
    for i=1:Nx
        if(x(i)<0)
            S(i)=Psi2(A1,T(1),x(i),k1);
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi(T(2),T(3),x(i),k2);
        elseif (x(i)>=b && x(i)<b+a)
            S(i)=Psi2(T(4),T(5),x(i),k3);
        elseif (x(i)>=b+a && x(i)<2*b+a)
            S(i)=Psi(T(6),T(7),x(i),k4);
        elseif (x(i)>=2*b+a)
            S(i)=Psi2(T(8),0,x(i),k5);
        end
    end
    plot(x*1e9,En/e+abs(S).^2/20,'b','LineWidth',3)
    hold on
%     plot(x*1e9,En/e+real(S)/20,'b--')
%     plot(x*1e9,En/e+imag(S)/20,'r-.')    
    xnew=[-b,   0,  0,  b,  b, a+b, a+b,    a+2*b, a+2*b, a+3*b];
    Unew=[0     0   U1  U1  0   0   U2      U2      0       0]/e;
    plot(xnew*1e9, Unew,'HandleVisibility','off')
    xlim([x(1) x(end)]*1e9)
    ylim([0 Umax])
    yline(En/e,'Color',[0.7 0.7 0.7]);
    xlabel('x,нм')
    ylabel('E, эВ')
%     legend('|\Psi|^2','Re(\Psi)','Im(\Psi)','E_0')
    grid on
    title('Динамика ВФ в профиле')

% 
%     if ((En==E0(1)) || (En==E0(2)))
%         Eline.Color=[0 0 0];
%         Eline.LineWidth=2;
%         writeVideo(v,getframe(gcf));
%         Eline.LineWidth=1;
%         writeVideo(v,getframe(gcf));
%         Eline.LineWidth=2;
%         writeVideo(v,getframe(gcf));
%     end

        
    m1=(m1+m2)/2; m2=m1;
    %solving Schr eq
    m=m1; mb=m2;
    gm1=@(E)sqrt(2*m*(E))/hbar;
    gm2=@(E)sqrt(2*mb*(-E+U1))/hbar;
    gm3=@(E)sqrt(2*m*(E))/hbar;
    gm4=@(E)sqrt(2*mb*(-E+U2))/hbar;
    gm5=@(E)sqrt(2*m*(E))/hbar;
    A1=1;
    
    MKoef=@(E)[ 1,              -1,                         -1,                                 0,                                  0,                                  0,                                   0,                                   0;
                -1i*gm1(E)/m,   -gm2(E)/mb,                 gm2(E)/mb,                          0,                                  0,                                  0,                                   0,                                   0;
                0,              exp(gm2(E)*b),              exp(-gm2(E)*b),                     -exp(1i*gm3(E)*b),                  -exp(-1i*gm3(E)*b),                 0,                                   0,                                   0;
                0,              gm2(E)/mb.*exp(gm2(E)*b),   -gm2(E)/mb.*exp(-gm2(E)*b),         -1i*gm3(E)/m.*exp(1i*gm3(E)*b),     1i*gm3(E)/m.*exp(-1i*gm3(E)*b),     0,                                   0,                                   0;
                0,              0,                          0,                                  exp(1i*gm3(E)*(a+b)),               exp(-1i*gm3(E)*(a+b)),              -exp(gm4(E)*(a+b)),                  -exp(-gm4(E)*(a+b)),                 0;
                0,              0,                          0,                                  1i*gm3(E)/m.*exp(1i*gm3(E)*(a+b)), -1i*gm3(E)/m.*exp(-1i*gm3(E)*(a+b)), -gm4(E)/mb.*exp(gm4(E)*(a+b)),       gm4(E)/mb.*exp(-gm4(E)*(a+b)),       0;
                0,              0,                          0,                                  0,                                  0,                                  exp(gm4(E)*(a+2*b)),                 exp(-gm4(E)*(a+2*b)),                -exp(1i*gm5(E)*(a+2*b));
                0,              0,                          0,                                  0,                                  0,                                  gm4(E)/mb.*exp(gm4(E)*(a+2*b)),      -gm4(E)/mb.*exp(-gm4(E)*(a+2*b)),    -1i*gm5(E)/m.*exp(1i*gm5(E)*(a+2*b))];
    MSvob=@(E)[-A1; -A1*1i*gm1(E)/m; 0; 0; 0; 0; 0; 0];
    % Koef=@(E)bicgstab(MKoef(E),MSvob(E), 1e-12, 100);
    % Koef=@(E)gmres(MKoef(E),MSvob(E),10,1e-12,15);
    Koef=@(E)(MKoef(E)\MSvob(E));
    
    Dav=zeros(1,N);
    for i=1:N
        Tav=Koef(E(i));
        Dav(i)=abs(Tav(8))^2;
    end
    Eav=E(islocalmax(Dav));
    Daver=Dav(islocalmax(Dav));
    
    Eaver(ind)=Eav(1);

     subplot(1,2,2)
    hold on
    plot(Dav,E/e,'--r')
    hold on
%     plot(D(id1),E(id1)/e,'--r','LineWidth',2)
%     plot(D(id2),E(id2)/e,'--b','LineWidth',2)
%     width=E(islocalmin(abs(D-1+exp(-1))));
%     width=[width(1),width(2)]/e;
%     plot((1-exp(-1))*[1,1],width,'k','LineWidth',2)
%     label = {sprintf(' E=%.0fмэВ \n D=%.1f \n g=%.1fмэВ',[E0(1)/e*1000,D0(1),(width(2)-width(1))*1000]),sprintf('E = %.0f мэВ\nD=%.1f',[E0(2)/e*1000,D0(2)])};
%     yl=yline(Eav(1)/e,'--r',label,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top');
    grid on
    ylabel('E,эВ')
    xlabel('D')
    title('Коэффициент прохождения')
    ylim([0 Umax])
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

    En=Eav(1);
    subplot(1,2,1)
    S=zeros(1,100);
    k1=gm1(En); k2=gm2(En); k3=gm3(En); k4=gm4(En); k5=gm5(En);
    Psi=@(A,B,x,k)A*exp(k*x)+B*exp(-k*x);
    Psi2=@(A,B,x,k)A*exp(1i*k*x)+B*exp(-1i*k*x);
%     N=300;
%     x=linspace(-b,a+3*b,N);
    T=Koef(En);
    % T=T./sqrt(T.*conj(T));
    for i=1:Nx
        if(x(i)<0)
            S(i)=Psi2(A1,T(1),x(i),k1);
        elseif (x(i)>=0 && x(i)<b)
            S(i)=Psi(T(2),T(3),x(i),k2);
        elseif (x(i)>=b && x(i)<b+a)
            S(i)=Psi2(T(4),T(5),x(i),k3);
        elseif (x(i)>=b+a && x(i)<2*b+a)
            S(i)=Psi(T(6),T(7),x(i),k4);
        elseif (x(i)>=2*b+a)
            S(i)=Psi2(T(8),0,x(i),k5);
        end
    end
    plot(x*1e9,En/e+abs(S).^2/20,'r','LineWidth',3)
    hold on
%     plot(x*1e9,En/e+real(S)/20,'b--')
%     plot(x*1e9,En/e+imag(S)/20,'r-.')    
%     xnew=[-b,   0,  0,  b,  b, a+b, a+b,    a+2*b, a+2*b, a+3*b];
%     Unew=[0     0   U1  U1  0   0   U2      U2      0       0]/e;
%     plot(xnew*1e9, Unew,'HandleVisibility','off')
%     xlim([x(1) x(end)]*1e9)
%     ylim([0 Umax])
    Eline=yline(En/e,'Color',[0.7 0.7 0.7]);
%     xlabel('x,нм')
%     ylabel('E, эВ')
%     legend('|\Psi|^2','Re(\Psi)','Im(\Psi)','E_0')
%     grid on
%     title('Динамика ВФ в профиле')


    writeVideo(v,getframe(gcf));
    close (f)
end
close (v)
close all

eps=abs(Eaver-Efct)./Efct;

%% graph
figure ('Units','normalized','OuterPosition', [0.05 0.05 0.45 0.45]);
plot(A,eps,'k')
grid on
ylabel('eps, %')
xlabel('a, нм')
title({'Погрешность экстремума в зависимости', 'от ширины канала','при ширине барьеров b=5 нм для Al_{0.3}Ga_{0.7}As'})

yyaxis right
plot(A,Efct/e,'-');
hold on
plot(A,Eaver/e,'--');
plot(A,abs(Eaver-Efct)/e,'-.')
legend('Относит.','Фактич. мас.', 'Средн. мас.', 'dE_{res}')
ylabel('E_{res}, эВ')
box on
xlim([A(1), A(end)])
  

s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
s.matlab.fonts.editor.normal.Size.TemporaryValue = 14;


% exportgraphics(f,'analyt.jpg')