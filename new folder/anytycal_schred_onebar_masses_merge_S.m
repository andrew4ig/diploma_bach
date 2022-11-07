clc
clear
close all

%constants
k=1.38e-23;
kT=k*300;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
eps0=8.85e-12;

zlim=[-6,10];


X=0.3;
%initial paramets
a=3e-9;                 %size of chanel, nm ?????????????
b=5e-9;                 %size of barrier, nm ?????????????
U0=e*X*0.74;            %height of barrier,  eV
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

Umax=1;
N=2500;
E=linspace(0,Umax*e,N);
Df=zeros(1,N);
for i=1:N
    T=Koef(E(i));
    Df(i)=abs(T(4))^2;
end

Emini=E(islocalmin(Df));
Emini=Emini(1);

Efact=E(E<Emini);
Df=Df(E<Emini);

%     if(isempty(Emaxi))
%         En=U0;
%     else
%         En=Emaxi(1);
%     end


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


Dav=zeros(1,N);
for i=1:N
    T=Koef(E(i));
    Dav(i)=abs(T(4))^2;
end

Emini=E(islocalmin(Dav));
Emini=Emini(1);

Eaver=E(E<Emini);
Dav=Dav(E<Emini);


figure ('Units','normalized','OuterPosition', [0.05 0.05 0.4 0.4]);
plot(Efact/e,Df,'*k','MarkerIndices',round(linspace(1,length(Df),80)))
hold on
plot(Eaver/e,Dav,'-.k')
ylabel('D')

yyaxis right
ax = gca;
ax.YColor = 'r';
if(length(Eaver)>length(Efact))
    En=Eaver;
else
    En=Efact;
end
mu=0.3;
V=0.05;
S=@(Ez)2*pi*m1*kT/(hbar^2)*log((1+exp((mu*e-Ez)/(kT)))./(1+exp((mu*e-Ez-V*e)/(kT))));
plot(En/e,S(En),'r');
hold on
plot(Efact/e,S(Efact).*Df, 'r*','MarkerIndices',round(linspace(1,length(Df),80)))
plot(Eaver/e,S(Eaver).*Dav,'r-.')
ylabel('S, ***')
xlabel('E, эВ')

legend('D_{facto}','D_{jure}', 'S','S*D_{facto}','S*D_{jure}','Location','best')
grid on
% title('Динамика ВФ в профиле')


fprintf('done')

s.matlab.fonts.editor.normal.Name.TemporaryValue = 'TimesNewRoman';
s.matlab.fonts.editor.normal.Size.TemporaryValue = 14;
