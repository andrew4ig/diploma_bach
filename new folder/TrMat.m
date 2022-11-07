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
a=3e-9;                 %size of chanel, 3nm
b=2e-9;                 %size of barrier, 2nm
U0=1*e;                 %height of barrier, 1 eV
m1=0.067*m0;             %eff. mass in GaAs, kg
m2=(0.067+0.083*0.3)*m0;%eff. mass in AlGaAs(0.3), kg
Ef=0.14*e;              %Fermi energy, (at T=300)
T=300;                  %room temperature, K
v=1e6;                  %frequency, 1Mhz
Ez=0.05*e;              %connection energy
V=@(t)sin(2*pi*v*t);    %periodic voltage, V


g=@(E,m,U)sqrt(2*m*(E-U))/hbar;
g1=@(E)g(E,m1,0); g2=@(E)g(E,m2,U0);
g=@(E)[g1(E),g2(E),g1(E),g2(E),g1(E)];
z=[0, b, b+a, 2*b+a];
m=[m1,m2,m1,m2,m1];

T1=@(E)[0.5*(1+g1(E)/g2(E)*m2/m1)*exp(-1i*(g2(E)-g1(E))*z(1)), 0.5*(1-g1(E)/g2(E)*m2/m1)*exp(-1i*(g2(E)+g1(E))*z(1));
       0.5*(1-g1(E)/g2(E)*m2/m1)*exp(1i*(g2(E)+g1(E))*z(1)), 0.5*(1+g1(E)/g2(E)*m2/m1)*exp(-1i*(g2(E)-g1(E))*z(1))];
T(1)
