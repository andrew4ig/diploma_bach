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

X=1;
%initial paramets
a=(3)*1e-9;                 %size of chanel, nm ?????????????
b=(2)*1e-9;                 %size of barrier, nm ?????????????
U0=X*0.74;               %height of barrier,  eV
m1=0.067*m0;                %eff. mass in GaAs, kg
m2=(0.067+0.083*X)*m0;      %eff. mass in AlGaAs(x), kg

dx1=1e-10;              %grid step
L1=2*b+a;
NE=1500;
N1=floor(L1/dx1);
x1=(0:N1-1)*dx1;          %x-vector

Vmax=1.5;
Volt=0:0.03:Vmax;
Umax=1.1*e;
Umax=max(Umax,U0*e);
En=linspace(0, Umax,NE); 

[J, D]=current(m1, m2, Volt, En, NE, a, b, U0, e, kT, hbar);
ma=(m1+m2)*0.5; mb=ma;
[Jav, Dav]=current(ma, mb, Volt, En, NE, a, b, U0, e, kT, hbar);

Jm=zeros(20,51);
Dm=zeros(20,NE);
eps=zeros(20,1);
index=0; 
%%
for i=0:0.05:1
    index=index+1
    ma=(m1*i+m2*(1-i));
    mb=ma;
    [Jm(index,:),Dm(index,:)]=current(ma, mb, Volt, En, NE, a, b, U0, e, kT, hbar);
    k_sim=(Jm(index,:)./J); 
    k_sim=mean(k_sim(~isnan(k_sim)));
    k_sim=(abs(J*k_sim-Jm(index,:))./J); 
    eps(index)=mean(k_sim(~isnan(k_sim)));
end

[z,y]=min(eps)

%%
v = VideoWriter('MassesAproach.mp4','MPEG-4');
v.FrameRate = 4;
v.Quality = 100;
open(v);
index=0;

for i=0:0.05:1
    index=index+1;
    ma=(m1*i+m2*(1-i));
    V=0;
    for ii=1:N1
            if(x1(ii)>=0&&x1(ii)<b)
    %                 mm(ii)=m1;
    %                 UU(ii)=-V*ii/N1;
                mm(ii)=m2;
                UU(ii)=U0-V*ii/N1;
            elseif (x1(ii)>=b&&x1(ii)<b+a)
                mm(ii)=m1;
                UU(ii)=-V*ii/N1;
    %                 mm(ii)=m2;
    %                 UU(ii)=U0-V*ii/N1;
            elseif(x1(ii)>=b+a)
    %                 mm(ii)=m1;
    %                 UU(ii)=-V*ii/N1;
                mm(ii)=m2;
                UU(ii)=U0-V*ii/N1;
            end
    end

    f=figure('Units','normalized','OuterPosition',[0.1 0.1 0.5 0.5]);
    subplot(1,3,1)
    x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
    y=[0 0 UU -V -V];
    Emax=En(islocalmax(D));
    Eavmax=En(islocalmax(Dm(index,:)));
    plot(x,y,'-b')
    x=x+dx1*1e9/2;
    hold on
%     plot(x,y,'_','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor',[0.95 0.95 0.95],'MarkerIndices',1:4:N1+4)
    if (~isempty(Emax))
    plot([b,a+b]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
    plot([b,a+b]*1e9,ones(2,length(Eavmax)).*Eavmax/e,'k','MarkerSize',2)
    end
    ylim([-Umax/e-.1, Umax/e])
    xlim([x(1)-0.1 x(end)+0.1])
    title({'???????????????????????????? ??????????????', '??????????????????'})
    xlabel('x,????')
    ylabel('E, ????')
    grid on

    yyaxis right
    plot(x1*1e9,mm/m0,'r')
    plot(x1*1e9,ma*ones(1,35)/m0,'r')
    text(3.5, 0.5, sprintf('m=m1*%.2f+\n+m2*(1-%.2f)', [i,i]),"FontSize",10,"FontWeight","bold",'HorizontalAlignment','center')
    ylabel('m/m_0')
    ylim([0 2])
    ax = gca;
    ax.YColor = 'r';
    set(gca,'FontSize',11,'fontWeight','bold')
    
    subplot(1,3,2)
    hold on
    fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
    plot(D,En/e,'k','LineWidth',1);
    plot(Dm(index,:),En/e,'--k','LineWidth',1);
    ylim([-Umax/e-.1, Umax/e])
    xlim([0 1])
    ylabel('E,????')
    xlabel('D')
    title({'??????????????????????' ,'????????????????????????'})
    grid on
    box on
    legend('????????. ??????????', "????. ??????????", 'Location','south')
    set(gca,'FontSize',11,'fontWeight','bold')
    
    
    subplot(1,3,3)
    hold on
    plot(Volt,J*1e6,'k','LineWidth',1);
    plot(Volt,Jm(index,:)*1e6,'--k','LineWidth',1);
    ylim([-2e-3, 5e-3])
%     ylim([0 5e-3])
    ylabel('J, ??????')
    xlabel('V, ??')
    title('??????')
    grid on
    box on
    legend('????????. ??????????', "????. ??????????", 'Location','south')
%         k_sim=(Jav./J); k_sim=mean(k_sim(~isnan(k_sim)));
%         k_sim=(abs(J*k_sim-Jav)./J); k_sim=mean(k_sim(~isnan(k_sim)));
    text((Volt(end)+Volt(1))/2, flip(ylim)*[0; 1]+(flip(ylim)*[1; -1]/4.5),...
        sprintf('eps = %0.2f', eps(index)),"FontSize",12,"FontWeight","bold",'HorizontalAlignment','center')
    set(gca,'FontSize',11,'fontWeight','bold')

    writeVideo(v,getframe(gcf));
    close(f)
end
close all
close(v)

%%
if isnan(J(1))
    J(1)=0;
    D(1)=0;
end
if isnan(Jav(1))
    Jav(1)=0;
    Dav(1)=0;
end

%% graphic
UU=zeros(1,N1);
V=0;
for ii=1:N1
   if(x1(ii)>=0&&x1(ii)<b)
%       UU(ii)=-V*ii/N1;
        UU(ii)=U0-V*ii/N1;
   elseif (x1(ii)>=b&&x1(ii)<b+a)
        UU(ii)=-V*ii/N1;
%       UU(ii)=U0-V*ii/N1;
   elseif(x1(ii)>=b+a)
%       UU(ii)=-V*ii/N1;
        UU(ii)=U0-V*ii/N1;
   end
end
UU=(UU-V/N1/2)*e;

f=figure('Units','normalized','OuterPosition',[0.1 0.1 0.45 0.45]);
subplot(1,3,1)
x=[-0.5 0 x1*1e9 x1(end)*1e9 x1(end)*1e9+0.5];
y=[0 0 UU/e -V -V];
Emax=En(islocalmax(D)); 
plot(x,y,'-b')
x=x+dx1*1e9/2;
hold on
plot(x,y,'_','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor',[0.95 0.95 0.95],'MarkerIndices',1:4:N1+4)
if (~isempty(Emax))
plot([b,a+b]*1e9,ones(2,length(Emax)).*Emax/e,'k--','MarkerSize',2)
end
% text((x1(end)+x1(1))/2, -Umax*0.9, sprintf('V = %0.2f ??', V),"FontSize",12,"FontWeight","bold",'HorizontalAlignment','center')
ylim([-Umax/e-.1, Umax/e])
xlim([x(1)-0.1 x(end)+0.1])
title({'???????????????????????????? ??????????????', '??????????????????'})
xlabel('x,????')
ylabel('E, ????')
grid on
set(gca,'FontSize',12,'fontWeight','bold')

subplot(1,3,2)
hold on
fill([0 0 1 1],[-Vmax-.1 0 0 -Vmax-.1],'r','FaceAlpha', 0.05,'HandleVisibility','off')
plot(D,En/e,'k','LineWidth',1);
plot(Dav,En/e,'--k','LineWidth',1);
ylim([-Umax/e-.1, Umax/e])
xlim([0 1])
ylabel('E,????')
xlabel('D')
title({'??????????????????????' ,'????????????????????????'})
grid on
box on
legend('????????. ??????????', "????. ??????????", 'Location','south')
%     k_sim=(Dav./D); k_sim=mean(k_sim(2:end));
%     disp=sqrt(sum(Dav-D*k_sim).^2/(length(NE)-1));
%     k_sim=(abs(D*k_sim-Dav)./D); k_sim=mean(k_sim(2:end));
% text(0.5, flip(ylim)*[0; 1]+(flip(ylim)*[1; -1]/4.5),...
%     sprintf('k_{??????????} = %0.2f', 1-k_sim),"FontSize",12,"FontWeight","bold",'HorizontalAlignment','center')
set(gca,'FontSize',12,'fontWeight','bold')


subplot(1,3,3)
hold on
plot(Volt,J*1e6,'k','LineWidth',1);
plot(Volt,Jav*1e6,'--k','LineWidth',1);
ylim([ylim*[1;0]-ylim*[-1; 1]/3, ylim*[0;1]])
% xlim([0 1])
ylabel('J, ??????')
xlabel('V, ??')
title('??????')
grid on
box on
legend('????????. ??????????', "????. ??????????", 'Location','south')
    k_sim=(Jav./J); k_sim=mean(k_sim(~isnan(k_sim)));
    k_sim=(abs(J*k_sim-Jav)./J); k_sim=mean(k_sim(~isnan(k_sim)));
text((Volt(end)+Volt(1))/2, flip(ylim)*[0; 1]+(flip(ylim)*[1; -1]/4.5),...
    sprintf('eps = %0.2f', k_sim),"FontSize",12,"FontWeight","bold",'HorizontalAlignment','center')
set(gca,'FontSize',12,'fontWeight','bold')

fprintf('done')

%% function calculating transparency and current
function [J,D1] = current(m1, m2, Volt, En, NE, a, b, U0, e, kT, hbar)
    dE=(En(2)-En(1));
    n=1e21;
    mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.3; 
    
    g=@(E,m,U)sqrt(2*m*(E-U))/hbar;
    T=@(m1,m2,g1,g2,z)[0.5*(1+g1/g2*m2/m1)*exp(-1i*(g2-g1)*z),   0.5*(1-g1/g2*m2/m1)*exp(-1i*(g2+g1)*z);
                       0.5*(1-g1/g2*m2/m1)*exp( 1i*(g2+g1)*z),   0.5*(1+g1/g2*m2/m1)*exp( 1i*(g2-g1)*z)];
    J=zeros(1,length(Volt));
%     D1=zeros(1,NE);

    for i=1:length(Volt)
        V=Volt(i);
        dx1=1e-10;              %grid step
        L1=2*b+a;               %total length
        N1=floor(L1/dx1);         %amount of grid cells
        x1=(0:N1-1)*dx1;          %x-vector
  
        mm=zeros(1,N1);
        UU=zeros(1,N1);
        
        for ii=1:N1
            if(x1(ii)>=0&&x1(ii)<b)
%                 mm(ii)=m1;
%                 UU(ii)=-V*ii/N1;
                mm(ii)=m2;
                UU(ii)=U0-V*ii/N1;
            elseif (x1(ii)>=b&&x1(ii)<b+a)
                mm(ii)=m1;
                UU(ii)=-V*ii/N1;
%                 mm(ii)=m2;
%                 UU(ii)=U0-V*ii/N1;
            elseif(x1(ii)>=b+a)
%                 mm(ii)=m1;
%                 UU(ii)=-V*ii/N1;
                mm(ii)=m2;
                UU(ii)=U0-V*ii/N1;
            end
        end
        
        UU=(UU-V/N1/2)*e;
        
        T0=@(E)T(m2,m1,g(E,m2,UU(N1)),g(E,m1,-V*e),x1(N1));
        for ii=N1-1:-1:1
            T0=@(E)T0(E)*T(mm(ii),mm(ii+1),g(E,mm(ii),UU(ii)),g(E,mm(ii+1),UU(ii+1)),x1(ii));
        end
        T0=@(E)T0(E)*T(m1,m2,g(E,m1,0),g(E,m2,UU(1)),x1(1)-x1(2));
        Tline=@(E)reshape(T0(E),1,[]);
    
        
        D=zeros(1,NE);
        for ii=1:NE
            T1=Tline(En(ii));
            D(ii)=abs(g(En,m1,-V*e))/abs(g(En,m1,0))*mm(1)/mm(5)*abs((T1(4)*T1(1)-T1(2)*T1(3))/T1(4))^2;
        end
       
        S=@(Ez)log((1+exp((mu*e-Ez)/(kT)))./(1+exp((mu*e-Ez-V*e)/(kT))));
        S=2*pi*m1*kT/(hbar^2)*S(En);
        J(i)=e/hbar*2/(2*pi)^3*trapz(S.*D)*dE*e; 
    
        if isnan(J(i)) && i>1
            J(i)=J(i-1)
        end

        if(i==1)
            D1=D;
        end
    end
end

