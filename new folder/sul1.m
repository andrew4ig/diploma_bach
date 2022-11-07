clear all 
NN = 50; 
hbar = 1.054e-34; 
m0 = 9.1e-31; 
melec = 0.25*m0; 
elec = 1.6e-19; 
eV2J = 1.6e-19;
J2eV = 1./eV2J; 
del_x = 2.e-10; % The cell size 
DX = del_x* 1e9;
XX = (DX:DX:NN* DX);
% Energies are eV 
chi0 = J2eV* hbar^2/(2*melec*del_x^2);
G0= elec^ 2/(2*pi*hbar); 
% Constant to calculate I 
kT = 0.025; 
% --- Channel potential 
V = zeros(1,NN);
%Resonant barrier 
for n= 9:13 %
    V(n) = 0.4; 
end
for n= 22:26 %
    V(n) = 0.4; 
end
for n= 35:39 %
    V(n) = 0.4; 
end
subplot(1,2,2) 
plot(XX,V) 
axis( [ 0 10 0 .5 ]) 
xlabel('nm') 
ylabel('V (eV)' ) 
set(gca,'fontsize',12)
Emax = 1; Emin = 0.; 
NE = 300; 
del_E = (Emax- Emin)/NE; 
EE = (0:del_E:del_E* (NE-1));
mu = 0.2;% input(' mu (left side) -->' ) 
% -- Calculate the Fermi function at the left contact ---
fermi1 =zeros(1,NE); 
for m= 1:NE 
    fermi1(m) = 1/(1 +  exp( (EE(m) -  mu ) /kT )); 
end
Npoints = 420; 
Iout = zeros(1,Npoints); 
Vin = zeros(1,Npoints);
 %----------- This is the main loop - ------------
for nvolt= 1:Npoints 
        VDS = 0.005* (nvolt-1)+ 0.0001; 
        Vin(nvolt) = VDS; 
        % --- Change in V for VDS - -----
        for n= 1:NN 
            VD(n) = -(n-1)*VDS/(NN-1); 
            V(n) = V(n) +  VD(n); 
        end
        Vmin = min(V);
    % -- Construct the Hamiltonian - -
    H = zeros(NN,NN); 
    H(1,1) = 2* chi0+V(1); 
    H(1,2) = -1*chi0; 
    for n= 2:NN-1 
        H(n,n-1)= -1*chi0; 
        H(n,n) = 2*chi0+ V(n); 
        H(n,n+1)= -1*chi0; 
    end
    H(NN,NN-1) = -1*chi0;
    H(NN,NN)   = 2*chi0+V(NN);
    % -- Calculate the Fermi function at the right contact - --
    fermi2 =zeros(1,NE); 
    TM =zeros(1,NE); 
    over =zeros(1,NE);
    for m= 1:NE 
        fermi2(m) = 1/(1 +  exp( (EE(m) -  (mu +  V(NN))) /kT )); 
    end
    % --- Calculate the Transmission function and the current 
    sigma1 = zeros(NN,NN); 
    sigma2 = zeros(NN,NN); 
    gamma1 = zeros(NN,NN); 
    gamma2 = zeros(NN,NN); 
    sig1 = 0.; 
    sig2 = 0.; 
    eta = 0; 
    n = zeros(NN,1); 
    I = 0.; 
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
        I = I +  (eV2J* del_E)*(elec/(2*pi*hbar))*TM(m)*(fermi1(m) - fermi2(m)); 
        %I = I +  G0* del_E*TM(m)*(fermi1(m) - fermi2(m)); 
        over(m) = TM(m)* (fermi1(m) - fermi2(m)); 
    end
    Iout(nvolt) = I;
end
Imax = 1e6* max(Iout); 
subplot(1,2,1) 
plot(Vin,1e6*Iout,'k') 
axis( [ 0 1 0 1.2* Imax ]) 
xlabel('V_D_S (eV)' ) 
ylabel('I (uA)' ) 
set(gca,'fontsize',12) 
title('Cal-IV') 
TT = text(.2,.7* Imax,'m_1'); 
set(TT,'fontsize',12) 
TT = text(.27,.72* Imax,sprintf(' = %4.2f eV' ,mu)); 
set(TT,'fontsize',12)
TT = text(.2,.4* Imax,'Resonant barrier' ); 
set(TT,'fontsize',12) 
grid on