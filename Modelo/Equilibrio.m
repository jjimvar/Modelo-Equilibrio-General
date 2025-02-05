function [f] = Equilibrio1(XX,PP)
%
y  = XX(1);
c  = XX(2);
in = XX(3);
k  = XX(4);
n  = XX(5);
w  = XX(6);
R  = XX(7);
%----------------------------------------------------------------
% 2. Parameters
%----------------------------------------------------------------
beta  = PP(1);
phi   = PP(2);
A     = PP(3);
theta = PP(4);
delta = PP(5);
g     = PP(6);
v     = PP(7);
gamma = PP(8);
tauC = PP(14);
tauN = PP(17);
tauK = PP(20);

%----------------------------------------------------------------
% 3. Model: Steady State Equations
%----------------------------------------------------------------
f(1) = -y + A*(k^theta)*(n^(1-theta));
f(2) = -y + c + in + g;
f(3) = -k + (1-delta)*k + in;
f(4) =  phi*(c^gamma)*n^(1/v) - w*(1-tauN)/(1+tauC); % Static equation KPR
f(5) =  w - (1-theta)*y/n;
f(6) =  R  - theta*y/k;
f(7) =  1 - beta*(1 - delta + R*(1-tauK));
f = f';
