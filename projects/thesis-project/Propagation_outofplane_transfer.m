%% Orbit Propagation + MPC control
% Esempio Out-of-plane transfer

clearvars
close all
clc

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Cose utili per la funzione JB2008:
global const PC % astronomical constants & planetary coefficients
global DTCdata eopdata %???
SAT_Const
constants
load DE440Coeff.mat
PC = DE440Coeff;

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for i=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

% read solar storm indices
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
global SOLdata
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read geomagnetic storm indices
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

addpath(genpath(fullfile(pwd, 'Utility Files')));

% Ts = 3;
mu = 3.986e5; % [km^3/s^2]
R_E = 6378; % [km] raggio medio della Terra
J_2 = 1.08263e-3; 

%% Condizioni di partenza esempio out-of-plane transfer tesi Belloni
e_c0=0.001;
a_c0= 6771; % [km]
p_c0 = a_c0*(1-e_c0^2); % semilatus rectum [km]

n_c0 = sqrt(mu/(a_c0)^3);  
T_orbit = 2*pi/n_c0;
Durata_sim = 7*T_orbit; % nella versione definitiva era stato messo 7

rho = 1.05259630795205e-3; % [kg/km^3]

i_c0 = 97.05 * pi/180; % inclinazione dell'orbita del chief
OM_c0 = 30 * pi/180;
om_c0 = 90 * pi/180;
M_c0        = 0 * pi/180; % Mean anomaly
nu_c0 = Tempo_inverso_aggiustato(M_c0, a_c0, e_c0);
r_c0 = p_c0/(1+e_c0*cos(nu_c0));
u_c0 = om_c0 + M_c0;

% Starting ROE in km
delta_a0 = 0;
delta_lambda0=0;
delta_e_x0=0;
delta_e_y0=0.2/a_c0;
delta_i_x0=0;
delta_i_y0=0.18/a_c0;

% Desired ROE in km
delta_aT = 0;
delta_lambdaT=0;
delta_e_xT=0;
delta_e_yT=0.2/a_c0;
delta_i_xT=0;
delta_i_yT=0.42/a_c0;

e_xc0 = e_c0*cos(om_c0);
e_yc0 = e_c0*sin(om_c0);
a_d0   = a_c0*(delta_a0+1);
e_xd0 = e_xc0+delta_e_x0;
e_yd0 = e_yc0+delta_e_y0;
e_d0 = sqrt(e_xd0^2+e_yd0^2);
om_d0  = atan2(e_yd0,e_xd0);
om_d0 = mod(om_d0, 2*pi);
i_d0 = i_c0 + delta_i_x0;
OM_d0 = OM_c0 + delta_i_y0/sin(i_c0);
u_d0 = u_c0 - delta_i_y0/tan(i_c0) + delta_lambda0;
M_d0 = u_d0 - om_d0;
nu_d0 = Tempo_inverso_aggiustato(M_d0, a_d0, e_d0);

OEc_0=[a_c0;e_c0;OM_c0;i_c0;om_c0;M_c0];
OE_d0=[a_d0;e_d0;OM_d0;i_d0;om_d0;M_d0];

mass_c=20;% [kg]
mass_d=20; % [kg]
C_Dc = 0; 
S_c = 0.1*(1e-3)^2; % [km^2]
B_c = C_Dc*S_c/mass_c; % [km^2/kg] coefficiente balistico del satellite chief
C_Dd = 2.1;
S_d = 0.1*(1e-3)^2; % [km^2]
B_d = C_Dd*S_d/mass_d; % [km^2/kg] coefficiente balistico del satellite deputy 1

C_R=1; % reflectance coefficient di chief e deputy

delta_B = (B_d-B_c)/B_d; % definizione usata nella tesi di Belloni (cioè
% quella che sto usando al momento)
% delta_B = B_d-B_c; % definizione di Guffanti (cioè quella che ho usato nella tesi)

ROE_T = [delta_aT,delta_lambdaT,delta_e_xT,delta_e_yT, delta_i_xT, delta_i_yT, delta_B]';

max_thrust=0.65*(1e-3)*(1e-3); % [kN]
min_thrust=0.35*(1e-3)*(1e-3); % [kN]
u_max = max_thrust/mass_d; % [km/s^2]
u_min = min_thrust/mass_d; % [km/s^2]

%% Simulazione usuale

clc
step=10; % time step usato nel file di confronto dei risultati
step_num=10; % time step usato nel file di propagazione numerica
Ts = 100;
% tic
fprintf('Run della simulazione iniziato in data: %s.\n', string(datetime('now')));
fprintf('Run della dinamica numerica in corso\n');
mdl = 'Dinamica_numerica_MPC.slx';
open_system(mdl)
tic; out=sim(mdl); toc;

%% Simulazione usuale: plot traiettoria numerica
num_x_rtn=zeros(length(out.tout),1);
num_y_rtn=zeros(length(out.tout),1);
num_z_rtn=zeros(length(out.tout),1);
num_x_rtn(:,1) = out.r_rtn.signals.values(1,1,:);
num_y_rtn(:,1) = out.r_rtn.signals.values(2,1,:);
num_z_rtn(:,1) = out.r_rtn.signals.values(3,1,:);

% inizio della parte aggiunta ora

num_u_R_ms2 = zeros(length(out.u_components_ms2.time),1);
num_u_T_ms2 = zeros(length(out.u_components_ms2.time),1);
num_u_N_ms2 = zeros(length(out.u_components_ms2.time),1);

num_u_R_ms2(:,1) = out.u_components_ms2.signals.values(:,1);
num_u_T_ms2(:,1) = out.u_components_ms2.signals.values(:,2);
num_u_N_ms2(:,1) = out.u_components_ms2.signals.values(:,3);

N_orbit_sample = out.u_components_ms2.time ./ T_orbit;

close all

% Plot componenti della spinta
% figure
% 
% c_green = [0.4660 0.6740 0.1880]; % verde MATLAB
% c_red   = [0.8500 0.3250 0.0980]; % rosso MATLAB
% c_blue  = [0.0000 0.4470 0.7410]; % blu MATLAB
% 
% stairs(N_orbit_sample, num_u_R_ms2.*1e3, 'LineWidth',2, 'Color',c_green)
% hold on
% stairs(N_orbit_sample, num_u_T_ms2.*1e3, 'LineWidth',2, 'Color',c_red)
% stairs(N_orbit_sample, num_u_N_ms2.*1e3, 'LineWidth',2, 'Color',c_blue)
% 
% yline( u_max.*1e6, '--k', 'u_{max}', 'LineWidth', 1.3)
% yline(-u_max.*1e6, '--k', '-u_{max}', 'LineWidth', 1.3)
% 
% hold off
% 
% grid on
% xlabel('Number of orbits [-]')
% ylabel('u [mm/s^2]')
% legend('u_R','u_T','u_N', 'Location','eastoutside')
% set(gcf,'Color','w');   % sfondo figura bianco
% set(gca,'Color','w');   % sfondo assi bianco

figure
set(gcf,'Color','w');

% ===== colori =====
c_green = [0.4660 0.6740 0.1880];
c_red   = [0.8500 0.3250 0.0980];
c_blue  = [0.0000 0.4470 0.7410];

% ===== plot =====
stairs(N_orbit_sample, num_u_R_ms2.*1e3, ...
       'LineWidth',1.3, 'Color',c_green)
hold on
stairs(N_orbit_sample, num_u_T_ms2.*1e3, ...
       'LineWidth',1.3, 'Color',c_red)
stairs(N_orbit_sample, num_u_N_ms2.*1e3, ...
       'LineWidth',1.3, 'Color',c_blue)

yline( u_max.*1e6, '--k', 'u_{max}',  'LineWidth',1.3)
yline(-u_max.*1e6, '--k', '-u_{max}', 'LineWidth',1.3)

hold off
grid on

xlabel('Number of orbits [-]')
ylabel('u [mm/s^2]')

legend('u_R','u_T','u_N','Location','eastoutside')

set(gca,'Color','w');

% ===== dimensione figura =====
set(gcf,'Units','centimeters');
set(gcf,'Position',[2 2 26 12]);   % larga e bassa

% ===== riduzione margini orizzontali =====
ax = gca;
ax.Units = 'normalized';
ax.Position = [0.08 0.15 0.78 0.78];
%              ↑      ↑
%     spazio grafico   lascia posto alla legenda a destra

% ===== impostazioni PDF =====
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 26 12]);
set(gcf,'PaperSize',[26 12]);

% ===== stampa finale =====
print(gcf,'control_profile_oopt.pdf','-dpdf','-painters');


%-------------------------------------------------------------

% Plot ellipsoidal altitude del satellite deputy

osc_ellis_h_km = zeros(length(out.tout),1);
osc_ellis_h_km(:,1) = out.osc_ellis_h_km.signals.values(:,1);
N_orbit = out.tout ./ T_orbit;

% figure
% set(gcf,'Units','centimeters');
% set(gcf,'Position',[2 2 18 12]);
% 
% plot(N_orbit, osc_ellis_h_km, 'k', 'LineWidth', 1.5)
% grid on
% 
% xlabel('Number of orbits [-]')
% ylabel('h [km]')
% 
% set(gcf,'Color','w');
% set(gca,'Color','w');
% 
% % ===== FIX REALE DEL TAGLIO YLABEL =====
% ax = gca;
% ax.Units = 'normalized';
% ax.Position = [0.16 0.12 0.83 0.80];
% %            ^^^^^
% %            margine sinistro più largo
% 
% exportgraphics(gcf,'ellipsoidal_altitude.pdf','ContentType','vector');

figure
set(gcf,'Color','w');

plot(N_orbit, osc_ellis_h_km, 'k', 'LineWidth', 1.5)
grid on

xlabel('Number of orbits [-]')
ylabel('h [km]')

set(gcf,'Units','centimeters');
set(gcf,'Position',[2 2 26 4]);

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 26 8]);
set(gcf,'PaperSize',[26 8]);

% ===== RIDUZIONE MARGINI ORIZZONTALI =====
ax = gca;
ax.Units = 'normalized';
ax.Position = [0.08 0.15 0.90 0.80];

print(gcf,'ellipsoidal_altitude_oopt.pdf','-dpdf','-painters');

% fine della parte aggiunta ora

n_d0 = sqrt(mu/(a_d0)^3);  
T_orbitd0 = 2*pi/n_d0;
a=a_d0; e=e_d0; incl=i_d0; OM=OM_d0; om=om_d0; nu=nu_d0;
mdl = 'RTNfreeorbit.slx';
open_system(mdl)
out2=sim(mdl);
x_rtn_startorb=zeros(length(out2.tout),1);
y_rtn_startorb=zeros(length(out2.tout),1);
z_rtn_startorb=zeros(length(out2.tout),1);
x_rtn_startorb(:,1) = out2.r_rtn.signals.values(1,1,:);
y_rtn_startorb(:,1) = out2.r_rtn.signals.values(2,1,:);
z_rtn_startorb(:,1) = out2.r_rtn.signals.values(3,1,:);

a_dT   = a_c0*(delta_aT+1);
e_xdT = e_xc0+delta_e_xT;
e_ydT = e_yc0+delta_e_yT;
e_dT = sqrt(e_xdT^2+e_ydT^2);
om_dT  = atan2(e_ydT,e_xdT);
om_dT = mod(om_dT, 2*pi);
i_dT = i_c0 + delta_i_xT;
OM_dT = OM_c0 + delta_i_yT/sin(i_c0);
u_dT = u_c0 - delta_i_yT/tan(i_c0) + delta_lambdaT;
M_dT = u_dT - om_dT;
nu_dT = Tempo_inverso_aggiustato(M_dT, a_dT, e_dT);

a=a_dT; e=e_dT; incl=i_dT; OM=OM_dT; om=om_dT; nu=nu_dT;
mdl = 'RTNfreeorbit.slx';
open_system(mdl)
out3=sim(mdl);
x_rtn_targetorb=zeros(length(out3.tout),1);
y_rtn_targetorb=zeros(length(out3.tout),1);
z_rtn_targetorb=zeros(length(out3.tout),1);
x_rtn_targetorb(:,1) = out3.r_rtn.signals.values(1,1,:);
y_rtn_targetorb(:,1) = out3.r_rtn.signals.values(2,1,:);
z_rtn_targetorb(:,1) = out3.r_rtn.signals.values(3,1,:);

% Plot posizioni relative nel frame RTN
% close all
figure 
subplot(1,3,1);  % 1 riga, 3 colonne, primo plot
plot(out.tout./(60*60), 1e3.*num_x_rtn);
% title('\underline{x} [m] vs ore simulazione');
title('$\underline{x}\ \mathrm{[m]}\ \mathrm{vs\ ore\ simulazione}$', ...
      'Interpreter','latex');
grid on
xlim([0 Durata_sim/(60*60)])
xlabel('$t\ \mathrm{[h]}$','Interpreter','latex');
ylabel('$\underline{x}\ \mathrm{[m]}$','Interpreter','latex');

subplot(1,3,2);  % secondo plot
plot(out.tout./(60*60), 1e3.*num_y_rtn);
title('y RTN numerica [m] vs ore simulazione');
grid on
xlim([0 Durata_sim/(60*60)])

subplot(1,3,3);  % terzo plot
plot(out.tout./(60*60), 1e3.*num_z_rtn);
title('z RTN numerica [m] vs ore simulazione');
grid on
xlim([0 Durata_sim/(60*60)])

set(gcf,'Color','w');   % sfondo figura bianco
set(gca,'Color','w');   % sfondo assi bianco

%% Plot orbita relativa tridimensionale in assi RTN

%-------------------------------
% codice chatgpt
figure
set(gcf,'Units','centimeters');
set(gcf,'Position',[2 2 20 15]);

plot3(1e3.*num_y_rtn, 1e3.*num_x_rtn, 1e3.*num_z_rtn, ...
      'k--', 'LineWidth', 0.9)   % <-- tratteggiata e più sottile

xlabel('$\underline{y}\ \mathrm{[m]}$','Interpreter','latex')
ylabel('$\underline{x}\ \mathrm{[m]}$','Interpreter','latex')
zlabel('$\underline{z}\ \mathrm{[m]}$','Interpreter','latex')

hold on
plot3(1e3.*num_y_rtn(1), 1e3.*num_x_rtn(1), 1e3.*num_z_rtn(1), ...
      'bs', 'MarkerSize', 9, 'MarkerFaceColor', 'b')
plot3(1e3.*num_y_rtn(end), 1e3.*num_x_rtn(end), 1e3.*num_z_rtn(end), ...
      'rs', 'MarkerSize', 9, 'MarkerFaceColor', 'r')

grid on
axis equal

plot3(1e3.*y_rtn_startorb, 1e3.*x_rtn_startorb, 1e3.*z_rtn_startorb, ...
      'b', 'LineWidth', 2)
plot3(1e3.*y_rtn_targetorb, 1e3.*x_rtn_targetorb, 1e3.*z_rtn_targetorb, ...
      'r', 'LineWidth', 2)

% title('Traiettoria satellite deputy in assi RTN');

legend('Deputy relative trajectory', 'Initial deputy position', 'Final deputy position', ...
       'Starting relative orbit', 'Target relative orbit');

set(gcf,'Color','w');
set(gca,'Color','w');

set(gca,'LooseInset', max(get(gca,'TightInset'), 0.04));
set(gcf,'Renderer','opengl');

exportgraphics(gcf,'3D_relative_trajectory_oopt.pdf','Resolution',600);

%-------------------------------------------------------------

% fine codice chatgpt

% punto iniziale della simulazione
p0 = [num_x_rtn(1), num_y_rtn(1), num_z_rtn(1)];

% orbita iniziale
orbit_start = [x_rtn_startorb(:), y_rtn_startorb(:), z_rtn_startorb(:)];

% distanza euclidea
distances_i = sqrt(sum((orbit_start - p0).^2, 2));

% distanza minima
min_dist_i = min(distances_i);

disp(['Distanza minima da punto iniziale simulazione a orbita teorica iniziale: ', num2str(min_dist_i), ' km'])

% punto finale della simulazione
p_end = [num_x_rtn(end), num_y_rtn(end), num_z_rtn(end)];
% matrice punti orbita target
orbit_target = [x_rtn_targetorb(:), y_rtn_targetorb(:), z_rtn_targetorb(:)];
% distanza euclidea
distances_t = sqrt(sum((orbit_target - p0).^2, 2));
% calcolo distanza minima
min_dist_t = min(distances_t);

disp(['Distanza minima da punto finale simulazione a orbita target: ', num2str(min_dist_t), ' km'])


%% Calcolo risultati per presentazione TASI

a_c_ROE_T_meters=struct;
a_c_ROE_T_meters.signals.dimensions=[6,1];
a_c_ROE_T_meters.time=out.tout;
a_c_osc_ROE_meters=struct;
a_c_osc_ROE_meters.signals.dimensions=[6,1];
a_c_osc_ROE_meters.time=out.tout;
u_components_ms2=struct;
u_components_ms2.signals.dimensions=[3,1];
u_components_ms2.time=out.tout;
osc_ellis_h_km=struct;
osc_ellis_h_km.signals.dimensions=[1,1];
osc_ellis_h_km.time=out.tout;

% Assegnazione valori

a_c_ROE_T_meters.signals.values(1:6,1,:)=out.a_c_ROE_T_meters.signals.values(:,1,:);
a_c_osc_ROE_meters.signals.values(1:6,1,:)=out.a_c_osc_ROE_meters.signals.values(:,1,:);
u_components_ms2.signals.values(:,1:3)=out.u_components_ms2.signals.values(:,1:3);
osc_ellis_h_km.signals.values(1,1,:)=out.osc_ellis_h_km.signals.values(1,1,:);

% Plot ROEs

tmax = max(out.tout(:));

% any(diff(out.tout) <= 0) % riga suggeritami da ChatGPT per verificare se il
% vettore che contiene i tempi è monotono. Se come risposta ottengo uno 0
% logico, va tutto bene; in caso invece ottenga un 1 logico, vuol dire che
% effettivamente c'è qualche elemento nel vettore che lo rende non
% monotono.

figure
set(gcf,'Color','w');

t = tiledlayout(3,2, ...
    'TileSpacing','loose', ...
    'Padding','compact');

for k = 1:6
    nexttile

    plot(N_orbit(:), ...
         squeeze(a_c_osc_ROE_meters.signals.values(k,1,:)), ...
         'k', 'LineWidth', 1.5);
    hold on
    yline(a_c_ROE_T_meters.signals.values(k,1), 'r', 'LineWidth', 1.5);
    hold off

    grid on
    xlim([0 (tmax+100)/T_orbit]);

    if k==1
        ylabel('a_{c0}\deltaa [m]');
    elseif k==2
        ylabel('a_{c0}\delta\lambda [m]');
    elseif k==3
        ylabel('a_{c0}\deltae_x [m]');
    elseif k==4
        ylabel('a_{c0}\deltae_y [m]');
    elseif k==5
        ylabel('a_{c0}\deltai_x [m]');
    elseif k==6
        ylabel('a_{c0}\deltai_y [m]');
    end

    xlabel('Number of orbits [-]');
end

lgd = legend({'Simulation value','Target value'}, ...
             'Orientation','horizontal');
lgd.Layout.Tile = 'north';

set(gcf,'Units','centimeters');
set(gcf,'Position',[2 2 26 22]);   % dimensione figura

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 26 22]);
set(gcf,'PaperSize',[26 22]);

print(gcf,'ROE_evolution_oopt.pdf','-dpdf','-painters');


%%


% num_mean_OE_c = struct;
% num_mean_OE_c.signals.dimensions = [6,1];
% num_mean_ROE  = struct;
% num_mean_ROE.signals.dimensions = [6,1];
% 
% % I seguenti vettori li chiamo "medi", ma saranno medi o osculanti a
% % seconda di come imposto gli switch nella simulazione della dinamica numerica
% num_mean_OE_c.signals.values = out.num_mean_OE_c.signals.values; 
% num_mean_ROE.signals.values = out.num_mean_ROE.signals.values;
% num_mean_OE_c.time = out.tout;
% num_mean_ROE.time = out.tout;
% 
% fprintf('Run della dinamica linearizzata in corso\n');
% mdl2 = 'Dinamica_linearizzata.slx';
% open_system(mdl2)
% out2=sim(mdl2);
% 
% fprintf('Run del confronto dei risultati in corso\n');
% mdl3 = 'Confronto_dinamiche.slx';
% open_system(mdl3)
% out3=sim(mdl3);
% 
% fprintf('Run della simulazione terminato in data: %s.\n', string(datetime('now')));
% durata_simulazione = toc;
% fprintf('Durata della simulazione: %d minuti.\n', durata_simulazione/60);


%% Simulazione con risultati integrazione numerica mediati alla Koenig

% clc
% step=10; % time step usato nel file di confronto dei risultati
% tic
% fprintf('Run della simulazione iniziato in data: %s.\n', string(datetime('now')));
% fprintf('Run della dinamica numerica in corso\n');
% 
% out = sim('Dinamica_numerica.slx');
% 
% passo = diff(out.tout);
% [~, q] = min(abs(out.tout - T_orbit));
% duration_m = out.tout(q);
% time_0 = out.tout;
% 
% num_mean_OE_c = struct;
% num_mean_OE_c.signals.dimensions = [6,1];
% num_mean_ROE  = struct;
% num_mean_ROE.signals.dimensions = [6,1];
% 
% num_mean_rho  = struct;
% num_mean_rho.signals.dimensions = 1;
% T_orbit_m=[];
% j=1;
% T_orbit_m(j)=T_orbit;
% T_rem=out.tout(end);
% i=1;
% 
% while (T_rem >= duration_m && duration_m >0)
% i_fin = q;
% a_mc=trapz(out.tout(i:i_fin), out.num_osc_OE_c.signals.values(1,1,i:i_fin)) ./ duration_m;
% e_mc=trapz(out.tout(i:i_fin), out.num_osc_OE_c.signals.values(2,1,i:i_fin)) ./ duration_m;
% i_mc=trapz(out.tout(i:i_fin), out.num_osc_OE_c.signals.values(3,1,i:i_fin)) ./ duration_m;
% OM_mc=trapz(out.tout(i:i_fin), out.num_osc_OE_c.signals.values(4,1,i:i_fin)) ./ duration_m;
% om_mc=trapz(out.tout(i:i_fin), out.num_osc_OE_c.signals.values(5,1,i:i_fin)) ./ duration_m;
% M_c=out.num_osc_OE_c.signals.values(6,:,i:i_fin);
% num_mean_OE_c.signals.values(1,:,i:i_fin) = a_mc;
% num_mean_OE_c.signals.values(2,:,i:i_fin) = e_mc;
% num_mean_OE_c.signals.values(3,:,i:i_fin) = i_mc;
% num_mean_OE_c.signals.values(4,:,i:i_fin) = OM_mc;
% num_mean_OE_c.signals.values(5,:,i:i_fin) = om_mc;
% num_mean_OE_c.signals.values(6,:,i:i_fin) = M_c;
% % Calcolo mean ROE
% mean_delta_a=trapz(out.tout(i:i_fin), out.num_osc_ROE.signals.values(1,1,i:i_fin)) ./ duration_m;
% mean_delta_lambda=trapz(out.tout(i:i_fin), out.num_osc_ROE.signals.values(2,1,i:i_fin)) ./ duration_m;
% mean_delta_e_x=trapz(out.tout(i:i_fin), out.num_osc_ROE.signals.values(3,1,i:i_fin)) ./ duration_m;
% mean_delta_e_y=trapz(out.tout(i:i_fin), out.num_osc_ROE.signals.values(4,1,i:i_fin)) ./ duration_m;
% mean_delta_i_x=trapz(out.tout(i:i_fin), out.num_osc_ROE.signals.values(5,1,i:i_fin)) ./ duration_m;
% mean_delta_i_y=trapz(out.tout(i:i_fin), out.num_osc_ROE.signals.values(6,1,i:i_fin)) ./ duration_m;
% num_mean_ROE.signals.values(1,:,i:i_fin) = mean_delta_a;
% num_mean_ROE.signals.values(2,:,i:i_fin) = mean_delta_lambda;
% num_mean_ROE.signals.values(3,:,i:i_fin) = mean_delta_e_x;
% num_mean_ROE.signals.values(4,:,i:i_fin) = mean_delta_e_y;
% num_mean_ROE.signals.values(5,:,i:i_fin) = mean_delta_i_x;
% num_mean_ROE.signals.values(6,:,i:i_fin) = mean_delta_i_y;
% % Calcolo mean rho
% m_rho = trapz(out.tout(i:i_fin), out.num_osc_rho.signals.values(i:i_fin,1)) ./ duration_m;
% num_mean_rho.signals.values(:,:,i:i_fin) = m_rho;
% 
% T_rem = Durata_sim - duration_m;
% i=i_fin;
% j=j+1;
% % T_orbit_m(j)=2*pi.*sqrt(out.num_osc_OE_c.signals.values(1,:,i).^3./mu);
% T_orbit_m(j) = T_orbit;
% [~, q] = min(abs(out.tout - j*T_orbit));
% duration_m = out.tout(q)-out.tout(i);
% 
% end
% 
% num_mean_OE_c.time = out.tout;
% num_mean_ROE.time = out.tout;
% 
% 
% fprintf('Run della dinamica linearizzata in corso\n');
% out2 = sim('Dinamica_linearizzata.slx');
% fprintf('Run del confronto dei risultati in corso\n');
% out3 = sim('Confronto_dinamiche.slx');
% fprintf('Run della simulazione terminato in data: %s.\n', string(datetime('now')));
% durata_simulazione = toc;
% fprintf('Durata della simulazione: %d minuti.\n', durata_simulazione/60);