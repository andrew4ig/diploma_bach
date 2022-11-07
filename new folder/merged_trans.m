%preparing workspace
clc
clear
close all

%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
kT=300*k;

X=0.3;
%initial paramets
U0=0.75;              %height of barrier,  eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*X)*m0*0+m1;  %eff. mass in AlGaAs(x), kg
m=m1;

V0=0.5;
V=linspace(0,V0,100);

%defining structure
j=2;                    %amount of barriers
dx=1e-10;               %grid step
a=3e-9;                 %size of chanel, nm
b=2e-9;                 %size of barrier, nm
L=j*b+(j+1)*a;          %total length
Np=floor(L/dx);         %amount of grid cells
koef=hbar^2/(2*m1*(dx^2))/e;
x=(1:Np)*dx;          %x-vector

n=1e21;
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.3; %chem potential

%correting hamiltonian corresponding to the heterstructure
E=eye(Np)*(2);
E=E+diag(ones(1,Np-1)*(-1),-1);
E=E+diag(ones(1,Np-1)*(-1),1);
E=E*koef;
for t=0:(2*j-1)
    if(mod(t,2)==0)
        E(:,round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L)+1)=E(:,round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L)+1)*(m1/m2);
    end
end
U=zeros(1,Np);
for t=0:(2*j-1)
    if(mod(t,2)==0)
        U(round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L)-1)=U0;
    end
end
H=E+diag(U);

One=[zeros(1,round(Np*a/L)-1), linspace(0,1,round(Np*(L-2*a)/L)+1),ones(1,round(Np*a/L))];   %empty vector

Nen=1000;
Umax=1.5;
Ener=linspace(0,Umax,Nen);

g=@(E,m,U)sqrt(2*m*(E-U))/hbar;
T=@(m1,m2,g1,g2,z)[0.5*(1+g1/g2*m2/m1)*exp(-1i*(g2-g1)*z),   0.5*(1-g1/g2*m2/m1)*exp(-1i*(g2+g1)*z);
                   0.5*(1-g1/g2*m2/m1)*exp( 1i*(g2+g1)*z),   0.5*(1+g1/g2*m2/m1)*exp( 1i*(g2-g1)*z)];

D=zeros(1,Nen);
D1=zeros(1,Nen);
J=zeros(1,100); JJ=J;

% v = VideoWriter('veryfVCC.avi','Motion JPEG AVI');
% v.Quality = 95;
% v.FrameRate=5;
% open(v);
isone=0;
ischecked=0;

for klo=1
    %correcting hamiltonian
    U1=-V(klo)*One;
    Ham=H+diag(U1);

    %Green's transmition coef-t
    S1=zeros(Np);S2=S1;
    for i=1:Nen
        En=Ener(i);
        S1(1,1)=-koef*exp(1i*sqrt(2*m*(En-U1(1))*e)/hbar*dx);
        S2(Np,Np)=-koef*exp(1i*sqrt(2*m*(En-U1(end))*e)/hbar*dx);

        G1=1i*(S1-S1');
        G2=1i*(S2-S2');

        Gr=(eye(Np)*(En)-Ham-S1-S2);

%         F=zeros(1,Np)'; A=53;
%         F(1)=2*1i*(-koef)*sin(sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
%         F(Np)=2*1i*(-koef)*sin(Np*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);
%         WaveFunction=Gr\F;
%         
        D(i)=real(trace(G1/Gr*G2/Gr'));
    end
   
    %transmission coefficient w/o votalge
    Emax=Ener(islocalmax(D));
    Dmax=D(islocalmax(D));
    if (isempty(D1))
        Emax=0;
        D1=0;
    end
    
    a1=3e-9;
    b1=2e-9;
    dx1=1e-10;              %grid step
    L1=2*b1+a1;               %total length
    N1=floor(L1/dx1);         %amount of grid cells
    x1=(0:N1-1)*dx1;          %x-vector
    % a=a+dx;b=b+dx;
    
    mm=zeros(1,N1);
    UU=zeros(1,N1);
    
    for ii=1:N1
    if(x1(ii)>=0&&x1(ii)<b1)
        mm(ii)=m2;
        UU(ii)=U0-V(klo)*ii/N1;
    elseif (x1(ii)>=b1&&x1(ii)<b1+a1)
        mm(ii)=m1;
        UU(ii)=-V(klo)*ii/N1;
    elseif(x1(ii)>=b1+a1)
        mm(ii)=m2;
        UU(ii)=U0-V(klo)*ii/N1;
    end
    end
    
    UU=(UU-V(klo)/N1/2)*e;
    
    T0=@(E)T(m2,m1,g(E,m2,UU(N1)),g(E,m1,-V(klo)*e),x1(N1));
    for ii=N1-1:-1:1
    T0=@(E)T0(E)*T(mm(ii),mm(ii+1),g(E,mm(ii),UU(ii)),g(E,mm(ii+1),UU(ii+1)),x1(ii));
    end
    T0=@(E)T0(E)*T(m1,m2,g(E,m1,0),g(E,m2,UU(1)),x1(1)-x1(2));
    Tline=@(E)reshape(T0(E),1,[]);
    Umax=1.5*e;%1*e*0.4;
    Umax=max(Umax,U0*e);    
    En=linspace(0,Umax,Nen);
    DD=zeros(1,Nen);
    for ii=1:Nen
    T1=Tline(En(ii));
    DD(ii)=abs(g(En,m1,-V(klo)*e))/abs(g(En,m1,0))*mm(1)/mm(5)*abs((T1(4)*T1(1)-T1(2)*T1(3))/T1(4))^2;
    end

    dE=En(2)-En(1);
    S=@(Ez)log((1+exp((mu*e-Ez)/(kT)))./(1+exp((mu*e-Ez-V(klo)*e)/(kT))));
    J(klo)=2*pi*m1*kT/(hbar^2)*e/hbar*2/(2*pi)^3*trapz(S(En).*D)*dE;
    JJ(klo)=2*pi*m1*kT/(hbar^2)*e/hbar*2/(2*pi)^3*trapz(S(En).*DD)*dE;

    DDmax=DD(islocalmax(DD));
   %%
    if(klo>0)
        f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.4 0.5]);
        set(gca,'FontSize',18)
        set(gca,'FontName','Times New Roman')
        
        if(ischecked==0&&length(Emax)==1)
            isone=1;
        end

        subplot(1,2,1)
        plot(x*1e9,U+U1)
        hold on
        ylim([-V0-0.1, 1.5])
        xlim([x(1), x(end)]*1e9)
        xlabel('x,нм')
        ylabel('U, эВ')
        grid on
        title('Энергетический профиль стуктуры')
        hold on
        if (~isempty(Emax))
        plot([a+b,2*a+b]*1e9,ones(2,length(Emax)).*Emax,'k--','MarkerSize',3)
        end
        text((x(end)+x(1))/2*1e9, -1, sprintf('V = %0.2f В', V(klo)),"FontSize",12,"FontWeight","bold",'HorizontalAlignment','center')


        load('analytival.mat')
        subplot(1,2,2)
        plot(D,Ener)
        ylabel('E,эВ')
        xlabel('D')
        title('Коэффициент прозрачности')
        ylim([-V0-0.1, 1.5])
        xlim([0 1])
        grid on
        hold on
        plot(Danalyt, Eanalyt/e,'r*','MarkerIndices',[1:50:300, 301:5:350,351:50:2000])
        plot(DD,Ener,'--k','MarkerSize',3)
        fill([0 0 1 1],[-V0-.1 0 0 -V0-.1],'r','FaceAlpha',0.05)
        legend('МФГ','ММП')

        if (1)
            dEner=0.1;
            rectangle('Position',[Dmax(1)*0.5 Emax(1)-dEner Dmax(1)*0.5 dEner*2])
            drawArrow = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'k','HandleVisibility','off' );
            x2 = [Dmax(1)/2 0.7];
            y2 = [Emax(1)-dEner -0.35];
            drawArrow(x2,y2)

            axes('Position',[0.74 0.2 0.15 0.14])
            box on
            indexs = Ener>Emax(1)-dEner & Ener<Emax(1)+dEner & (D > Dmax(1)*0.5 | DD > DDmax(1)*0.5);
            hold on
            plot(D(indexs),Ener(indexs))
            plot(DD(indexs),Ener(indexs),'k--','MarkerSize',2)
            plot(Danalyt(Eanalyt/e<Emax(1)+dEner &Eanalyt/e>Emax(1)-dEner), Eanalyt(Eanalyt/e<Emax(1)+dEner &Eanalyt/e>Emax(1)-dEner)/e,'r*')
            ylim([min(Ener(indexs)) max(Ener(indexs))])
            xlim([min(Dmax(1), DDmax(1))*0.5, max(DDmax(1),Dmax(1))])
            buf=ylim;
%             axis tight
%             %ylim([Emax(1)-0.05 Emax(1)+0.05])
            buf=(ceil(10000*linspace(buf(1),buf(2),5))/10000);
            buf(1)=ylim*[1; 0]; buf(end)=ylim*[0; 1];
            yticks(buf);
            buf=xlim;
            buf=(ceil(10000*linspace(buf(1),buf(2),5))/10000);
            buf(1)=xlim*[1; 0]; buf(end)=xlim*[0; 1];
            xticks(buf);
            grid minor
            
        end
%         writeVideo(v,getframe(gcf));
%         close(f)
    end
end
% close(v)

%%
f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.4 0.5]);
plot(V, J/1)
hold on
plot(V, JJ,'k--')
grid on
legend('МФГ','ММП')
xlabel('V, В')
ylabel('I, А/м^2')
grid on
title('ВАХ')
[~,b]=max(abs(JJ-J));
xline(V(b),'--')
plot([V(b),V(b)])


% %%
% E1=Ener(islocalmax(DD)); E1=E1(1);%расположение первого РУ по ММП
% 
% E2=Ener(islocalmax(D)); E2=E2(1); %расположение первого РУ по МФГ
% 
% E3=Eanalyt(islocalmax(Danalyt))/e; E3=E3(1); %расположение первого РУ по аналитическому методу
% 
% (E1-E3)/E1*100                                                                                   %#ok<NOPTS> 
% (E2-E3)/E2*100                                                                                                                  %#ok<NOPTS> 