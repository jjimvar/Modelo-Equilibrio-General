% Basic RBC Model 
%
close all;
%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------
var y c g in k n w R z tauC tauN tauK tau;
varexo epsilon omega omega_c omega_n omega_k;
parameters beta phi A theta delta G v rho_z sigma_z rho_f sigma_f psi, gamma;
parameters tauC_SS rho_c sigma_c tauN_SS rho_n sigma_n tauK_SS rho_k sigma_k;

%----------------------------------------------------------------
% 2. Calibration (from mcer12.m)
%----------------------------------------------------------------
PP      = load('parameters.txt');
beta    = PP(1);
phi     = PP(2);
A       = PP(3);
theta   = PP(4);
delta   = PP(5);
G       = PP(6);
v       = PP(7);
gamma   = PP(8);
rho_z   = PP(9);
sigma_z = PP(10);
rho_f   = PP(11);
sigma_f = PP(12);
psi     = PP(13);
tauC_SS = PP(14);
rho_c    = PP(15);
sigma_c  = PP(16);
tauN_SS = PP(17);
rho_n    = PP(18);
sigma_n  = PP(19);
tauK_SS = PP(20);
rho_k    = PP(21);
sigma_k  = PP(22);
tau_SS   = PP(23);

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model; 
  y   = A*exp(z)*(k(-1)^theta)*(n^(1-theta));
  c*(1+tauC) + in = w*n*(1-tauN) + R*k(-1)*(1-tauK) - tau;
  g   = (1-rho_f)*G + rho_f*g(-1) + psi*(y-1) + sigma_f*omega;
  g   = c*tauC + tauN*w*n + tauK*R*k(-1) + tau;
  k   = (1-delta)*k(-1) + in;
  phi*(c^gamma)*n^(1/v) = w*(1-tauN)/(1+tauC);
  w   = (1-theta)*y/(n);
  R   = theta*y/k(-1);
  c^(-gamma) = beta*((1+tauC)/(1+tauC(+1)))*(1 - delta + R(+1)*(1-tauK(+1)))*c(+1)^(-gamma);
  z    = rho_z*z(-1) + sigma_z*epsilon;
  tauC = (1-rho_c)*tauC_SS + rho_c*tauC(-1) + omega_c*sigma_c;
  tauN = (1-rho_n)*tauN_SS + rho_n*tauN(-1) + omega_n*sigma_n;
  tauK = (1-rho_k)*tauK_SS + rho_k*tauK(-1) + omega_k*sigma_k;
end;

%----------------------------------------------------------------
% 4. Computation (from mcer12.m)
%----------------------------------------------------------------
YY = load('steady_state.txt');

initval;
  y  = YY(1);
  c  = YY(2);
  in = YY(3);
  g  = YY(4);
  k  = YY(5);
  n  = YY(6);
  w  = YY(7);
  R  = YY(8);
  z  = 0; 
  epsilon = 0;
  omega   = 0;
  omega_c = 0;
  omega_n = 0;
  omega_k = 0;
  tau     = tau_SS;
  tauC    = tauC_SS;
  tauN    = tauN_SS;
  tauK    = tauK_SS;
end;

shocks;
  var epsilon = 1;
  var omega   = 1;
  var omega_c = 1;
  var omega_n = 1;
  var omega_k = 1;
end;

steady;
%----------------------------------------------------------------
% 5. Some Results
%----------------------------------------------------------------
stoch_simul(periods=5000,irf=40,order=1,drop=200);