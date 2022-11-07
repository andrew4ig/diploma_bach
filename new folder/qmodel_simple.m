function qmodel_simple()

% система исток - уровень в яме - сток

qe = 1.60E-19;
k=1.38e-23;
T=300;
hbar=1.0546e-34;
m = 0.063 * 9.11E-31;
e0 = 0.3 * qe;
l = 2E-9;
f = 1E12;

hs = @(x)x>=0;
df = @(x)log(1+exp((0.12*qe-x)/0.03/qe)).*hs(x);
es = @(n)3E-1*qe*n;
e = @(v) e0-qe*v/2;
v = @(t) 1/2*sin(2*pi*f*t-pi/2)+1/2;

    function [dydt, j] = rtdsystem(v, n)
        er = e(v)+es(n);
        pw = sqrt(2*m*e0)/m;
        ps = sqrt(2*m*(er))/m;
        pd = sqrt(2*m*(er+qe*v))/m;
        wsw = ps./l .* hs(er);
        wws = pw./l .* hs(er);
        wdw = pd./l .* hs(er+qe*v);
        wwd = pw./l .* hs(er+qe*v);
        osw = wsw .* (1 - n) .* df(er);
        odw = wdw .* (1 - n) .* df(er+qe*v);
        ows = wws .* n;
        owd = wwd .* n;
        dydt = (osw + odw) - (ows + owd);
        j = (osw + owd) - (odw + ows);
        j = qe*j;
    end

    function [n, j, T] = nd(v, t) 
            function dydt = odefun(t, n)
                dydt = rtdsystem(v(t), n);
            end
        [T, n] = ode15s(@odefun, t, 0, odeset('RelTol',1E-10));
        [~, j] = rtdsystem(v(t), n);
    end

% t = linspace(0,10/f,10000)';
% [n, j] = nd(@(t)+v(t), t);
% figure;
% plot(t, n);

t = linspace(0,2/f,10000)';
[n, j, T] = nd(@(t)+v(t), t);

figure;
subplot(1,4,1)
plot(v(T), j);
ylabel('$J,{A \over m^2*s}$','Interpreter','latex')
xlabel('$Voltage,V$','Interpreter','latex')
title('$VCC$','Interpreter','latex')

subplot(1,4,2)
plot(T, n);
xlabel('$time,ms$','Interpreter','latex')
ylabel('$n$','Interpreter','latex')
title('$n=n(t)$','Interpreter','latex')

subplot(1,4,3)
plot(T, j);
xlabel('$time,ms$','Interpreter','latex')
ylabel('$J,{A \over m^2*s}$','Interpreter','latex')
title('$J=J(t)$','Interpreter','latex')

subplot(1,4,4)
plot(T, v(T));
xlabel('$time,ms$','Interpreter','latex')
ylabel('$Voltage,V$','Interpreter','latex')
title('$V=V(t)$','Interpreter','latex')

% vt = v(t(9000:end));
% j = j(9000:end);
% vp = vt(1:floor(end/2));
% jp = j(1:floor(end/2));
% vm = vt(floor(end/2):end);
% jm = j(floor(end/2):end);
% 
% figure;
% plot(vp, jp, vm, jm);
% axis([0 1 0 1E-4]);
% xlabel('Напряжение, В');
% ylabel('Ток');



end