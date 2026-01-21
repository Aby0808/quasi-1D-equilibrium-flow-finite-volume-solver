clear
clc
close all

% main script
% quasi 1D nozzle flow using influence coefficients


x0 = 0;
xt = 0.0499;
xe = 0.225;

r0 = 0.05;
rt = 0.0109;
re = 0.0432;

P0 = 42249;
T0 = 1000;

rold = r0;
Aold = pi*rold^2;
Mold = 0.01;
Pold = P0;
Told = T0;
N = 400;

Aall = zeros(N,1);
Mall = zeros(N,1);
Pall = zeros(N,1);
Tall = zeros(N,1);

Aall(1) = Aold;
Mall(1) = Mold;
Pall(1) = Pold;
Tall(1) = Told;

c = 1;

for x = (xe-x0)/N:(xe-x0)/N:xe
    if x < xt
        r = r0 + (rt-r0)*x/(xt-x0);
    elseif x > xt
        r = rt + (re-rt)*(x-xt)/(xe-xt);
    end
    A = pi*r^2;

    dA = A-Aold;
    [gamma, ~, ~] = FUN_get_gamma(Told, Pold);

    dM2 = -2*(1 + 0.5*(gamma-1)*Mold^2)*(dA*Mold^2 /Aold)/(1-Mold^2);
    dP = (gamma*Mold^2)*(dA*Pold /Aold)/(1-Mold^2);
    dT = ((gamma-1)*Mold^2)*(dA*Told /Aold)/(1-Mold^2);

    M = sqrt(Mold^2 + dM2);
    P = Pold + dP;
    T = Told + dT;

    Mold = M;
    Pold = P;
    Told = T;

    c=c+1;
    Aall(c) = A;
    Mall(c) = M;
    Pall(c) = P;
    Tall(c) = T;
end

plot(x,Aall);
xlabel('x(m)')
ylabel('Area (m^2)')
grid on

plot(x,Mall);
xlabel('x(m)')
ylabel('Mach')
grid on

plot(x,Pall);
xlabel('x(m)')
ylabel('Pressure (Pa)')
grid on

plot(x,Tall);
xlabel('x(m)')
ylabel('Temperature (K)')
grid on