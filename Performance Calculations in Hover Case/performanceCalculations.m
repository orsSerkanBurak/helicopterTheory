clc; clear; close all;
g = 9.81; % Gravitional Acceleation [m/sec^2]
rho_SL = 1.225; % Density of air at sea level [kg/m^3]	
alt = [1000,3000,5000];% Operating density altitude values [m]
P_SL = 2*530*10^3; % Engine power at sea level [W]
R = 11.93/2;% Main rotor radius [m]
S = pi*R^2;% Actuator disk area [m^2]
K = 1.15;% Induced power coefficient
W = 4000*9.81;% Weight of helicopter [N]
T = W;% Thrust of main rotor in hover case [N]
omega = 300*(1/60)*(2*pi);% Main Rotor Angular Speed [rad/sec]
Cd0 = 0.0012;% Drag coefficient
c = 0.405;% Chord of main rotor blades [m]
Nb = 4;% Number of main rotor blades
solidity = Nb*c/(pi*R);% Solidity of main rotor
X_tr = 11.63; % Distance of tail rotor from center of the main rotor [m]
R_tr = 1.1/2;% Radius of tail rotor blades [m]
S_tr = pi*R_tr^2;% Tail rotor disk area [m^2]
omega_tr = 1500*(1/60)*(2*pi); % Angular velocity of tail rotor [rad/sec]
c_tr = 0.2;% Chord of tail rotor blades [m]
Nb_tr = 13;% Number of tail rotor blades
solidity_t = Nb_tr*c_tr/(pi*R_tr);% Solidity of tail rotor blades
%% Part (a) 
A = alt*((0.00357)/(0.3048*518.4));
rho = ((1-A).^(1/0.235))*rho_SL;
%% Part (b)
P = P_SL.*((1-A).^(1/0.235));
%% Part (c)
P_main_i = (K.*T.^1.5)./(sqrt(2.*rho*S));
P_main_0 = S.*rho.*((omega.*R)^3).*solidity.*Cd0/8;
P_main_h = P_main_i + P_main_0;
C_p_main = P_main_h./(rho.*S.*((omega*R)^3));
%% Part (d)
Q_main_h = P_main_h./omega;
T_tr = Q_main_h./X_tr;
%% Part (e)
P_tail_i = (K*T_tr.^1.5)./(sqrt(2.*rho*S_tr));
P_tail_0 = S_tr.*rho.*((omega_tr*R_tr)^3)*solidity_t*Cd0/8;
P_tail_h = P_tail_i + P_tail_0 ;
C_p_tr = P_tail_h./(rho.*S_tr.*((omega_tr*R_tr)^3));
%% Part (f)
delta_p = P - (P_main_h+P_tail_h); % Excess power at each altitude
%% Part (g): 
C_t =T./(rho*S*(omega*R)^2);
Cpi = C_t.^(3/2)./sqrt(2);
FM  = Cpi./(K*Cpi + (solidity*Cd0/8)); % Calculating the Figure of Merit (FM) at each altitude
%% Part (h)
fig = uifigure;
table = uitable(fig);
d = {'Air density [kg/m^3]',rho(1),rho(2),rho(3);
'Available power [kW]',P(1)*10^-3,P(2)*10^-3,P(3)*10^-3;
'Required induced power of main rotor [kW]',P_main_i(1)*10^-3,P_main_i(2)*10^-3,P_main_i(3)*10^-3;
'Required profile power of main rotor [kW]',P_main_0(1)*10^-3,P_main_0(2)*10^-3,P_main_0(3)*10^-3;
'Required induced power of tail rotor [kW]',P_tail_i(1)*10^-3,P_tail_i(2)*10^-3,P_tail_i(3)*10^-3;
'Required profile power of tail rotor [kW]',P_tail_0(1)*10^-3,P_tail_0(2)*10^-3,P_tail_0(3)*10^-3;
'Excess power [kW]',delta_p(1)*10^-3,delta_p(2)*10^-3,delta_p(3)*10^-3;
'Figure of merit',FM(1),FM(2),FM(3);};
table.Data = d;
table.ColumnName = {'Altitude','1 km','3 km','5 km'};
table.Position = [20 20 520 400];
