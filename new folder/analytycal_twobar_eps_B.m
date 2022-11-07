clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;

B=0.5:0.2:10;
Ef=zeros(1,length(B));
Ea=zeros(1,length(B));

for ind=1:length(B)
    b=B(ind)*1e-9;
    b2=b;

    X=0.3;
    X2=0.3;
    %initial paramets
    a=3e-9;                 %size of chanel, nm ?????????????
%     b=2e-9;                 %size of barrier, nm ?????????????
%     b2=b;
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
    
    Umax=1.5;%max(U1,U2)*2/e;
    N=2000;
    E=linspace(0,Umax*e,N);
    D=zeros(1,N);
    for i=1:N
        T=Koef(E(i));
        D(i)=abs(T(8))^2;
    end
    Emax=E(islocalmax(D));
    Dmax=D(islocalmax(D));
    
    Eupper=linspace(Emax(1)*0.95,Emax(1)*1.05,N/2);
    Dupper=zeros(1,N/2);
    for i=1:N/2
        T=Koef(Eupper(i));
        Dupper(i)=abs(T(8))^2;
    end
    
    En=Eupper(islocalmax(Dupper));

    Ef(ind)=En(1);
   
    m=(m1+m2)/2; mb=m;
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
    
    Eav=linspace(0,Umax*e,N);
    Dav=zeros(1,N);
    for i=1:N
        T=Koef(Eav(i));
        Dav(i)=abs(T(8))^2;
    end
    Emaxav=Eav(islocalmax(Dav));
    Dmaxav=Dav(islocalmax(Dav));
    
    Eupperav=linspace(Emaxav(1)*0.95,Emaxav(1)*1.05,N/2);
    Dupperav=zeros(1,N/2);
    for i=1:N/2
        T=Koef(Eupperav(i));
        Dupperav(i)=abs(T(8))^2;
    end
    
    Enav=Eupperav(islocalmax(Dupperav));

    Ea(ind)=Enav(1);
end

eps=abs(Ef-Ea)./Ef;

%% graph
figure ('Units','normalized','OuterPosition', [0.05 0.05 0.45 0.45]);
plot(B*1e9,eps,'k')
grid on
ylabel('eps, %')
xlabel('b, нм')
title({'Погрешность экстремума в зависимости', 'от ширины барьеров'})

yyaxis right
plot(B*1e9,Ef/e,'-');
hold on
plot(B*1e9,Ea/e,'--');
plot(B*1e9,abs(Ea-Ef)/e,'-.')
legend('Относит.','Фактич. мас.', 'Средн. мас.', 'dE_{res}')
ylabel('E_{res}, эВ')
box on
xlim([B(1), B(end)]*1e9)


s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
s.matlab.fonts.editor.normal.Size.TemporaryValue = 14;