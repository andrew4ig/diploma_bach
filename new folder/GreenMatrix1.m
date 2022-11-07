%preparing workspace
clc
clear
close all


%constants
k=1.38e-23;
hbar=1.0546e-34;
m0=9.1e-31;
e=1.6e-19;
T=300;

X=1;
%initial paramets
a=3e-9;                 %size of chanel, nm
b=2e-9;                 %size of barrier, nm
U0=X*0.74;              %height of barrier,  eV
m1=0.067*m0;            %eff. mass in GaAs, kg
m2=(0.067+0.083*X)*m0;  %eff. mass in AlGaAs(x), kg

Nen=3000;
V0=1;
V=linspace(0,V0,Nen);

%defining structure
j=2;                    %amount of barriers
dx=0.4e-10;               %grid step
% a=a+1*dx;b=b+6*dx;
L=j*b+(j+1)*a;          %total length
Np=floor(L/dx);         %amount of grid cells
koef=hbar^2/(2*m1*(dx^2))/e;
x=(1:Np)*dx;          %x-vector

Jwf=@(F)(F.'*gradient(conj(F))-F'*gradient(F));
% Jwf=@(F)(F*gradient(F')-conj(F)*gradient(F.'));

n=1e21;
mu=hbar^2/(2*m1)*(3*pi^2*n)^(2/3)/e*0+0.3; %chem potential
% mu=mu*(1-pi^2/12*(k*300/e/mu)^2)

%correting hamiltonian corresponding to the heterstructure
E=eye(Np)*(2);
E=E+diag(ones(1,Np-1)*(-1),-1);
E=E+diag(ones(1,Np-1)*(-1),1);
E=E*koef;
for t=0:(2*j-1)
    if(mod(t,2)==0)
        E(round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L),:)=E(round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L),:)/(m1/m2);
    end
end
U=zeros(1,Np);
for t=0:(2*j-1)
    if(mod(t,2)==0)
        U(round(Np*(floor(t/2)*(b+a)+a)/L):round(Np*(floor(t/2)+1)*(b+a)/L))=U0;
    end
end
H=E+diag(U);

% b=b+6*dx;
One=[zeros(1,floor(Np*b/L)), linspace(0,1,floor(Np*(L-2*b)/L+1e-5)),ones(1,floor(Np*b/L))];   %empty vector
% b=b-6*dx;

Umax=1;
Ener=linspace(0,Umax,Nen);

D=zeros(1,Nen);
D1=zeros(1,Nen);
J=zeros(1,length(V));
JJ=zeros(1,length(V));

% v = VideoWriter('Greencorr.avi','Motion JPEG AVI');
% v.Quality = 95;
% v.FrameRate=5;
% open(v);

WFMAX=zeros(2,Np);
test=[200:200:Nen];
testid=0;

klo=1;

%correcting hamiltonian
U1=-V(klo)*One;
Ham=H+diag(U1);

%Green's transmition coef-t
S1=zeros(Np);S2=S1;
for i=1:Nen
    En=Ener(i);
    S1(1,1)=-koef*exp(1i*sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
    S2(Np,Np)=-koef*exp(1i*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);

    G1=1i*(S1-S1');
    G2=1i*(S2-S2');

    Gr=(eye(Np)*(En)-Ham-S1-S2);
    F=zeros(1,Np)'; A=1;
    F(1)=2*1i*(-koef)*sin(sqrt(2*m1*(En-U1(1))*e)/hbar*dx);
    F(Np)=2*1i*(-koef)*sin(Np*sqrt(2*m1*(En-U1(end))*e)/hbar*dx);
    WaveFunction=Gr\F;
    WaveFunction=WaveFunction*sqrt(WaveFunction.'*conj(WaveFunction));

    D(i)=real(trace(G1/Gr*G2/Gr'));

    if sum(i==test)

        testid=testid+1;
        WFMAX(testid,:)=WaveFunction';

        WF1=WaveFunction(x<a);
        WFLast=WaveFunction(x>2*b+2*a);
%             figure('Units','normalized','OuterPosition', [0.1 0.1 0.3 0.5]);
%             plot(x',abs(WaveFunction))
%             hold on
%             plot(x(x<a)',abs(WF1),'--')
%             plot(x(x>2*b+2*a)',abs(WFLast),'--')
%             hold off
%             J1(i)=Jwf(WF1);
%             JL(i)=Jwf(WFLast);
        D1(i)=abs(Jwf(WFLast))/A^2;%/abs(Jwf(WF1));
        
        RE=real(WF1);
        IM=imag(WF1);
        I=abs(WF1);
        X=x(x<a)';

        k=sqrt(2*m1*(En-U1(1))*e)/hbar;
        
        fitresult = cell( 2, 1 );
        gof = struct( 'sse', cell( 2, 1 ), ...
            'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
        
        
        [xData, yData] = prepareCurveData( X, IM );
        
%             Set up fittype and options.
        ft = fittype( 'b0+b1*sin(w*x)+b2*sin(2*w*x)+b3*sin(3*w*x)+b4*sin(4*w*x)+b5*sin(5*w*x)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf k];
        opts.StartPoint = [0.505957051665142 0.699076722656686 0.890903252535798 0.959291425205444 0.547215529963803 0.196595250431208  k];
        opts.Upper = [Inf Inf Inf Inf Inf Inf k];
        
        % Fit model to data.
        [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
        
        
%         [xData, yData] = prepareCurveData( X, RE );
        
%             Set up fittype and options.
        ft = fittype( 'a0+a1*cos(w*x)+a2*cos(2*w*x)+a3*cos(3*w*x)+a4*cos(4*w*x)+a5*cos(5*w*x)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf  k];
        opts.StartPoint = [0.149294005559057 0.257508254123736 0.840717255983663 0.254282178971531 0.814284826068816 0.929263623187228  k];
        opts.Upper = [Inf Inf Inf Inf Inf Inf  k];
        
        % Fit model to data.
        [fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );
        
        
        aa=[fitresult{2}.a0, fitresult{2}.a1, fitresult{2}.a2, fitresult{2}.a3, fitresult{2}.a4, fitresult{2}.a5];
        bb=[fitresult{1}.b0, fitresult{1}.b1, fitresult{1}.b2, fitresult{1}.b3, fitresult{1}.b4, fitresult{1}.b5];

        A=(aa+bb)/2;
%             A(1)=[];
%             sum(A.^2);
%             Jwf(WFLast)
        Ddiv(testid)=abs(i*hbar/2/m1*Jwf(WFLast))/(sum(A.^2));
        Ddiv1(testid)=abs(Jwf(WFLast))/(sum(A.^2));
        Ddiv2(testid)=abs(i*hbar/2/m1*Jwf(WFLast))/(sum(A.^2))/((WF1.'*conj(WF1)));
        Entest(testid)=En;
    end
end
if(isnan(D1(1)))
    D1(1)=0;
end

%%transmission coefficient w/o votalge
% if(klo==1)
%     Emax=Ener(islocalmax(D));
%     D1=D;
%     if (isempty(D1))
%         Emax=0;
%         D1=0;
%     end
% end


f=figure ('Units','normalized','OuterPosition', [0.05 0.05 0.5 0.5]);

subplot(1,2,1)
plot(x*1e9,U, 'k')
hold on
grid on
ylabel('E,эВ')
xlabel('x, нм')
title('Динамика ВФ в профиле')
norm=0.05*(WFMAX(1,:)*WFMAX(1,:)');
plot(x*1e9,abs(WFMAX(1,:)).^2/norm+Ener(test(1)),'-r','LineWidth',2)
plot(x*1e9,real(WFMAX(1,:))/sqrt(norm)*0.2+Ener(test(1)),'--r','MarkerIndices',round(linspace(1,length(x),50)))
plot(x*1e9,imag(WFMAX(1,:))/sqrt(norm)*0.2+Ener(test(1)),'.r','MarkerIndices',round(linspace(1,length(x),50)))
yline(Ener(test(1)),'Color',[0.7 0.7 0.7])
norm=0.05*(WFMAX(2,:)*WFMAX(2,:)');
plot(x*1e9,abs(WFMAX(2,:)).^2/norm+Ener(test(2)),'-b','LineWidth',2)
plot(x*1e9,real(WFMAX(2,:))/sqrt(norm)*0.2+Ener(test(2)),'--b','MarkerIndices',round(linspace(1,length(x),50)))
plot(x*1e9,imag(WFMAX(2,:))/sqrt(norm)*0.2+Ener(test(2)),'.b','MarkerIndices',round(linspace(1,length(x),50)))
yline(Ener(test(2)),'Color',[0.7 0.7 0.7])
xlim([0 x(end)]*1e9)
ylim([0 Umax])
hold off

subplot(1,2,2)
plot(D,Ener)
hold on; 
grid on
ylabel('E,эВ')
xlabel('D')
title('Коэффициент туннелирования')
xlim([0 1])
