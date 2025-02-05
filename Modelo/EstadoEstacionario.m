clear all
close all
%----------------------------------------------------------------
% In this exercise, utility can logarithmic or KPR:
%    u(c,1-N) = ln(c) + phi·ln(1-N)               ==> log-log case
%    u(c,1-N) = ln(c) - phi·(N^(1+1/v))/(1+1/v).  ==> KPR case
% and the production function is Cobb-Douglas:
%             y = A·(k^theta)·(N^(1-theta)).
%----------------------------------------------------------------
% Targeted moments:
y       = 1.00;      % Steady state output
g       = 0.25218;      % Steady state government spending
in      = 0.21311;      % Steady state investment
n       = 0.30;      % Steady-State fraction of hours worked
k       = (3.71905*4)*y; % Steady state capital/output
c       = y-g-in;    % Steady state consumption
% Exogenous parameters:
gamma   = 2;          % CRRA coefficient
theta   = 0.34512;    % Capital income share
v       = 0.72*2;     % Frisch elasticity
rho_z   = 0.84;       % Persistency of productivty shocks
sigma_z = 0.007;      % Standard Dev. of productivty shocks
psi     = -0.01/1.1;  % Automatic stabilizer

% Public Spending parameters
rho_f   = 0.9479;     % Persistency of fiscal shocks
sigma_f = 0.0003*3.5; % Standard Dev. of fiscal shocks
tauC_SS = 0.18948;
rho_c   = 0.810628;
sigma_c = 0.006617/4;
tauN_SS = 0.386869;
rho_n   = 0.854631;
sigma_n = 0.017335/4;
tauK_SS = 0.347096;
rho_k   = 0.793012;
sigma_k = 0.02776/4;

%----------------------------------------------------------------
% Parameters that require solving the model:
w     = (1-theta)*y/n;                        % Real wage
phi   = (w*(1-tauN_SS)/(1+tauC_SS))/((c^gamma)*(n^(1/v))); % MRS leisure-c when utility is KPR
delta = in/k;                                 % Steady state capital accumulation
R     = theta*y/k;                            % Demand for capital
r     = -1 + (R*(1-tauK_SS) + 1 - delta);       % Real interest rate
beta  = 1/(1+r);                              % Subjective discount rate
A     = y/((k^theta)*n^(1-theta));            % Technological scale parameter
G     = g;                                    % Steady state Gov. spending
tau   = G - tauC_SS*c - tauN_SS*w*n - tauK_SS*R*k;
PP    = [beta, phi, A, theta, delta, G, v, gamma, rho_z, sigma_z,...
    rho_f, sigma_f, psi, tauC_SS, rho_c, sigma_c, tauN_SS, rho_n, sigma_n,...
    tauK_SS, rho_k, sigma_k, tau]; % Save parameters
save parameters.txt PP -ascii -double
YY = [y,c,in,g,k,n,w,R]; % Save Steady-State values
save steady_state.txt YY -ascii -double
%----------------------------------------------------------------
% Check and Model sensitivity:
XX   = [y,c,in,k,n,w,R]'; %Steady state
para = PP;
f    = Equilibrio1(XX,PP); % Comprobamos que XX es un equilibrio estacionario
% para(1) = beta*0.99;  % Disminuye el descuento subjetivo;
% para(2) = phi*1.1; % Aumentamos permanente la desutilidad del trabajo;
% para(3) = A*1.1; % Aumento permanente 10% PTF;
% para(4) = theta + 0.07; % Aumemto de la fracción de rentas de capital;
% para(5) = delta*2;    % Duplicamos la tasa de depreciación del capital
% para(6) = G + 0.03; % Aumentamos el gasto G = 0.17 + 0.03;
% para(7) = v*2; % Duplicamos la elasticidad Frisch;
xx0 = XX; % Semilla == Seed
xx  = NEWTON_Equilibrio1(xx0,para);
Cambio = [xx0, xx, xx/xx(1)];
Nombres = {'PIB, y','Consumo, c','Inversión, in','Capital, k',...
    'Horas trabajadas, n', 'Salario real, w', 'Precio de arrendamiento, R'}';
Casos = {'Referencia','Nuevo','Relativo_al_PIB'};
Tab1  = array2table(round(Cambio,3),'RowNames',Nombres,'VariableNames',Casos);
Tab1
%----------------------------------------------------------------
% Execute dynare
dynare rbc1
% % REMEMBER: 
% % 1. Cada serie es simulada bajo el nombre definido en rbc1.mod 
% %    (véase la linea 7: var y c g in k n w R z f;)
% % 2. Las funciones de impulso respuesta (FIR) están recogidas con el nombre 
% %    de la variable y el choque fundamental. Por ejemplo, la FIR del
% %    consumo frente a un choque fiscal, se guarda en c_omega; la FIR de la
% %    inversión frente a un choque de productividad en in_epsilon; etc.
% % 3. Nada sale bien desde el principio. Todo es cuestión de esfuerzo.
% %    Ánimo.
% %----------------------------------------------------------------
% Express simulated variables in % fluctuations:
yy = log(y/mean(y)); % 
cc = log(c/mean(c)); % When c/0.63 = (1+zc) (zc == % dev. from Css=0.63), log(1+gc) ~ gc.
ii = log(in/mean(in));
gg = log(g/mean(g));
nn = log(n/mean(n));
ww = log(w/mean(w));
RR = log(R/mean(R));
Cycle = [yy, cc, ii, gg, nn];
Cy    = Cycle(:,1);
[T,k] = size(Cycle);
Tabla = [];
ret   = 4;
ade   = 4;
for j = 1:k
    BC    = Cycle(:,j);
    BC1   = correl(Cy,BC,ret,ade);
    Tabla = [Tabla; BC1'];
end
Lag_4  = Tabla(:,1);
Lag_3  = Tabla(:,2);
Lag_2  = Tabla(:,3);
Lag_1  = Tabla(:,4);
Center = Tabla(:,5);
Lead_1 = Tabla(:,6);
Lead_2 = Tabla(:,7);
Lead_3 = Tabla(:,8);
Lead_4 = Tabla(:,9);
Sigma_Abs = std(Cycle)';
Sigma_Rel = Sigma_Abs/Sigma_Abs(1);
% Correlograma = table(Etiquetas,Sigma_Abs, Sigma_Rel,...
%                      Lag_4,Lag_3,Lag_2,Lag_1,Center,...
%                      Lead_1,Lead_2,Lead_3,Lead_4);
Etiquetas = {'Pib', 'Consumo', 'Inversion', 'Gasto Publico', 'Empleo'};
Tabla = [Sigma_Abs, Sigma_Rel, Tabla];
CorreSim = array2table(round(Tabla,3),'RowNames',Etiquetas);
CorreSim
% % %----------------------------------------------------------------
% % % Simulated correlogram
% QQ = [NaN NaN -3:3];
% C1  = correl(yy,yy,3,3);
% C2  = correl(yy,cc,3,3);
% C3  = correl(yy,ii,3,3);
% C4  = correl(yy,gg,3,3);
% C5  = correl(yy,nn,3,3);
% C6  = correl(yy,ww,3,3);
% C7  = correl(yy,RR,3,3);
% 
% TABLA1 = [QQ; std(yy) std(yy)/std(yy) C1'; 
%               std(cc) std(cc)/std(yy) C2'; 
%               std(ii) std(ii)/std(yy) C3'; 
%               std(gg) std(gg)/std(yy) C4';
%               std(nn) std(nn)/std(yy) C5'; 
%               std(ww) std(ww)/std(yy) C6'; 
%               std(RR) std(RR)/std(yy) C7'];
% Nombres2 = {'a ','PIB, y','Consumo, c','Inversión, in','Gasto público, g',...
%     'Horas trabajadas, n','Salario real, w', 'Precio de arrendamiento, R'}';
% Codigos = {'SD','SD_rel','Lag3','Lag2','Lag1','Center','Lead1','Lead2','Lead3'};
% Tab2    = array2table(round(TABLA1,3),'RowNames',Nombres2,'VariableNames',Codigos);
% Tab2
% %----------------------------------------------------------------
% % Variance decomposition:
% Vy  = [cumsum(y_epsilon.^2)./cumsum(y_epsilon.^2 + y_omega.^2), cumsum(y_omega.^2)./cumsum(y_epsilon.^2 + y_omega.^2)];
% Vc  = [cumsum(c_epsilon.^2)./cumsum(c_epsilon.^2 + c_omega.^2), cumsum(c_omega.^2)./cumsum(c_epsilon.^2 + c_omega.^2)];
% Vin = [cumsum(in_epsilon.^2)./cumsum(in_epsilon.^2 + in_omega.^2), cumsum(in_omega.^2)./cumsum(in_epsilon.^2 + in_omega.^2)];
% Vg  = [cumsum(g_epsilon.^2)./cumsum(g_epsilon.^2 + g_omega.^2), cumsum(g_omega.^2)./cumsum(g_epsilon.^2 + g_omega.^2)];
% Vn  = [cumsum(n_epsilon.^2)./cumsum(n_epsilon.^2 + n_omega.^2), cumsum(n_omega.^2)./cumsum(n_epsilon.^2 + n_omega.^2)];
% Vw  = [cumsum(w_epsilon.^2)./cumsum(w_epsilon.^2 + w_omega.^2), cumsum(w_omega.^2)./cumsum(w_epsilon.^2 + w_omega.^2)];
% %----------------------------------------------------------------
% % Expenditure multipliers:
% % YY = [y,c,in,g,k,n,w,R];
% yss      = YY(1); % Steady state values
% css      = YY(2); % Steady state values
% iss      = YY(3); % Steady state values
% gss      = YY(4); % Steady state values
% kss      = YY(5); % Steady state values
% nss      = YY(6); % Steady state values
% wss      = YY(7); % Steady state values
% Rss      = YY(8); % Steady state values
% T        = length(g_omega);
% rr       = (1/beta + R_omega).^((0:T-1)');
% yy_omega = y_omega/yss;
% cc_omega = c_omega/css;
% ii_omega = in_omega/iss;
% gg_omega = g_omega/gss;
% kk_omega = k_omega/kss;
% nn_omega = n_omega/nss;
% ww_omega = w_omega/wss;
% RR_omega = R_omega/Rss;
% zz_omega = z_omega;
% yy_epsilon = y_epsilon/yss;
% cc_epsilon = c_epsilon/css;
% ii_epsilon = in_epsilon/iss;
% gg_epsilon = g_epsilon/gss;
% kk_epsilon = k_epsilon/kss;
% nn_epsilon = n_epsilon/nss;
% ww_epsilon = w_epsilon/wss;
% RR_epsilon = R_epsilon/Rss;
% zz_epsilon = z_epsilon;
% 
% PVy = ((cumsum(yy_omega./rr))./cumsum(gg_omega./rr))*(yss/gss);
% PVc = ((cumsum(cc_omega./rr))./cumsum(gg_omega./rr))*(css/gss);
% PVi = ((cumsum(ii_omega./rr))./cumsum(gg_omega./rr))*(iss/gss);
% PVg = ((cumsum(gg_omega./rr))./cumsum(gg_omega./rr))*(gss/gss);
% PV  = [PVy, PVc, PVi, PVg];
% PV  = [[1 4 8 T]', PV([1 4 8 end],:)];
% Nombres3 = {'Impact','1 year','2 year','Long term'}';
% Codigos  = {'Horizon','GDP','Conumption','Investment','Expenditure'};
% Tab3 = array2table(round(PV,3),'RowNames',Nombres3,'VariableNames',Codigos);
% Tab3
% %
% close all
% figure(1)
% subplot(3,3,1)
% plot(yy_epsilon,'linewidth',2)
% grid
% title('GDP IRF, $$\epsilon \rightarrow GDP$$','Interpreter','latex')
% 
% subplot(3,3,2)
% plot(cc_epsilon,'linewidth',2)
% grid
% title('Consumption IRF, $$\epsilon \rightarrow C$$','Interpreter','latex')
% 
% subplot(3,3,3)
% plot(ii_epsilon,'linewidth',2)
% grid
% title('Investment IRF, $$\epsilon \rightarrow I$$','Interpreter','latex')
% 
% subplot(3,3,4)
% plot(gg_epsilon,'linewidth',2)
% grid
% axis([0 40 -0.01 0.01])
% title('Expenditure IRF, $$\epsilon \rightarrow G$$','Interpreter','latex')
% 
% subplot(3,3,5)
% plot(nn_epsilon,'linewidth',2)
% grid
% title('Labor IRF, $$\epsilon \rightarrow N$$','Interpreter','latex')
% 
% subplot(3,3,6)
% plot(kk_epsilon,'linewidth',2)
% grid
% title('Capital IRF, $$\epsilon \rightarrow K$$','Interpreter','latex')
% 
% subplot(3,3,7)
% plot(ww_epsilon,'linewidth',2)
% grid
% title('Wage IRF, $$\epsilon \rightarrow w$$','Interpreter','latex')
% 
% subplot(3,3,8)
% plot(RR_epsilon,'linewidth',2)
% grid
% title('Interest rate IRF, $$\epsilon \rightarrow r$$','Interpreter','latex')
% 
% subplot(3,3,9)
% plot(zz_epsilon,'linewidth',2)
% grid
% title('TFP IRF, $$\epsilon \rightarrow z$$','Interpreter','latex')
% 
% figure(2)
% subplot(3,3,1)
% plot(yy_omega,'linewidth',2)
% grid
% title('GDP IRF, $$\omega \rightarrow GDP$$','Interpreter','latex')
% 
% subplot(3,3,2)
% plot(cc_omega,'linewidth',2)
% grid
% title('Consumption IRF, $$\omega \rightarrow C$$','Interpreter','latex')
% 
% subplot(3,3,3)
% plot(ii_omega,'linewidth',2)
% grid
% title('Investment IRF, $$\omega \rightarrow I$$','Interpreter','latex')
% 
% subplot(3,3,4)
% plot(gg_omega,'linewidth',2)
% grid
% title('Expenditure IRF, $$\omega \rightarrow G$$','Interpreter','latex')
% 
% subplot(3,3,5)
% plot(nn_omega,'linewidth',2)
% grid
% title('Labor IRF, $$\omega \rightarrow N$$','Interpreter','latex')
% 
% subplot(3,3,6)
% plot(kk_omega,'linewidth',2)
% grid
% title('Capital IRF, $$\omega \rightarrow K$$','Interpreter','latex')
% 
% subplot(3,3,7)
% plot(ww_omega,'linewidth',2)
% grid
% title('Wage IRF, $$\omega \rightarrow w$$','Interpreter','latex')
% 
% subplot(3,3,8)
% plot(RR_omega,'linewidth',2)
% grid
% title('Interest rate IRF, $$\omega \rightarrow r$$','Interpreter','latex')
% 
% subplot(3,3,9)
% plot(zz_omega,'linewidth',2)
% grid
% title('TFP IRF, $$\omega \rightarrow z$$','Interpreter','latex')
% %
% tiempo = toc/60