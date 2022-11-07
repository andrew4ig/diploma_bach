clear
figure('Units','normalized','OuterPosition',[0 0 1 1])
NN = 100; 
hbar = 1.06e-34; 
m0 = 9.11e-31; 
melec = 0.063*m0; 
eccoul = 1.6e-19;
eV2J = 1.6e-19; 
J2eV = 1./eV2J;

del_x = 1.e-10; 
DX = del_x* 1e9; 
X = (DX:DX:NN* DX);

N_center = NN/2; % Energies are eV 
chi0 = J2eV* hbar^2/(2*melec*del_x^2) ;

% ---- Specify the potential ----
V = zeros(1,NN); 

for n= 16:35 % 
    V(n) = 1; 

end
for n= 66:85 % 
    V(n) = 1; 
end

subplot(1,2,1) 
plot(X,V,'k') 
title('Trans') 
% axis( [ 0 10 -.1 .6 ]) 
grid on 
xlabel(' x (nm)' ) 
ylabel('V (eV)' ) 
ylim([0 1.5])


H = zeros(NN,NN); 
H(1,1) = 2* chi0+V(1);
H(1,2) = -1*chi0; 
for n= 2:NN-1 
    H(n,n-1)= -1*chi0; 
    H(n,n) = 2*chi0+ V(n); 
    H(n,n+1)= -1*chi0; 
end
H(NN,NN-1) = -1*chi0; 
H(NN,NN) = 2*chi0+V(NN); 

% -- Specify the energy range ---
Emax = 1.5; 
Emin = 0.; 
NE = 1000; 

del_E = (Emax- Emin)/NE; 
EE = (0:del_E:del_E* (NE-1));
sigma1 = zeros(NN,NN); 
sigma2 = zeros(NN,NN); 
gamma1 = zeros(NN,NN); 
gamma2 = zeros(NN,NN); 
sig1 = 0.; 
sig2 = 0.; 

eta = 1e-12; 
TM = zeros(1,NN); 

for m= 1:NE 
    k = sqrt(2* melec*(EE(m)-V(1))*eV2J)/hbar; 
    sig1 = exp(1i* k*del_x); 
    k = sqrt(2* melec*(EE(m)-V(NN))*eV2J)/hbar; 
    sig2 = exp(1i* k*del_x); 
    sigma1(1,1) = -chi0*sig1; 
    sigma2(NN,NN) = -chi0*sig2; 
    gamma1 = 1i* (sigma1-sigma1'); 
    gamma2 = 1i* (sigma2-sigma2'); 
    G = inv( (EE(m) + 1i* eta)*eye(NN) - H -  sigma1 -  sigma2); 
    TM(m) = real(trace(gamma1* G*gamma2*G')); 
end

subplot(1,2,2) 
plot(TM,EE,'k') 
% Elmax=EE(islocalmax(TM));
% yline(Elmax(1),'--',['first peak ', num2str((Elmax(1)*1000)), ' meV']);
grid on 
ylabel('E (eV)' ) 
xlabel('TM')
ylim([Emin Emax])