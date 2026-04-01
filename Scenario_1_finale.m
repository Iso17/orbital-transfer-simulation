close all
clc
clear


load("DatiSC1-2025.mat");
th0=0;
thn=2*pi;
dth= 2*pi/200;

%Parametri cartesiani orbita finale 
r_f=[-1524.866
     9996.648	
     4556.442];        % [km]
v_f=[-6.421	
     -2.655
     3.6760];          % [km/s]

%Parametri orbitali orbita iniziale
mu=398600;              % [km^3/s^2]
a=24400;                % [km]
e=0.7283;               % [-]
i=0.1047;               % [rad]
OM=0.8671;              % [rad]
om=3.107;               % [rad]
th_i=3.098;             % [rad]

%Parametri cartesiani orbita iniziale
[r_i, v_i] = par2car (a, e, i, OM, om, th_i, mu);
r_i_mod = norm(r_i);
v_i_mod = norm(v_i);

% Parametri orbitali orbita finale
[a_f,e_f,i_f,OM_f,om_f,th_f] = car2par(r_f, v_f, mu);

%% MANOVRA STANDARD PA

%---------------- manovra cambio di piano (deltat1)
[DeltaV_COP_PA, omf_PA, theta_PA] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);
deltat_COP_PA = TOF (a, e, th_i, theta_PA, mu);

th=th0:dth:thn;
N_PA=length(th);
rr_P_PA=zeros(3,N_PA);
vv_P_PA=zeros(3,N_PA);
for j=1:length(th)
[rr_P_PA(:,j),vv_P_PA(:,j)]=par2car(a,e,i,OM,om,th(j),mu);
end

rr_CO_PA=zeros(3,N_PA);
vv_CO_PA=zeros(3,N_PA);
for j=1:length(th)
[rr_CO_PA(:,j),vv_CO_PA(:,j)]=par2car(a,e,i_f,OM_f,omf_PA,th(j),mu);
end

dth1_PA= (theta_PA+2*pi)/200;
th1_PA=th_i:dth1_PA:(theta_PA+2*pi);
N1_PA=length(th1_PA);
rr_P1_PA=zeros(3,N1_PA);
vv_P1_PA=zeros(3,N1_PA);
for j=1:length(th1_PA)
[rr_P1_PA(:,j),vv_P1_PA(:,j)]=par2car(a,e,i,OM,om,th1_PA(j),mu);
end

[r_COPi_PA,~] = par2car(a, e, i, OM, om, th_i, mu); %punto di inizio
[r_COP_PA,~] = par2car(a, e, i, OM, om, theta_PA, mu); %punto di manovra

figure(Name='manovra cambio di piano') 
hold on
plot3(rr_P_PA(1,:),rr_P_PA(2,:),rr_P_PA(3,:))
plot3(rr_CO_PA(1,:),rr_CO_PA(2,:),rr_CO_PA(3,:))
plot3(r_COPi_PA(1), r_COPi_PA(2), r_COPi_PA(3), 'o', 'MarkerSize',6,'Color',	[0, 0.4470, 0.7410]);
plot3(r_COP_PA(1), r_COP_PA(2), r_COP_PA(3), 'o', 'MarkerSize',6,'MarkerFaceColor',[0, 0.4470, 0.7410]);
plot3(rr_P1_PA(1,:),rr_P1_PA(2,:),rr_P1_PA(3,:), 'LineWidth',2,'Color',[0, 0.4470, 0.7410])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
legend('Orbita iniziale', 'Orbita cambio di piano','punto iniziale','punto di manovra COP')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)

%--------------------- cambio pericentro (deltat2) 

% Prima intersezione 
[DeltaV_CPA_PA, thi_PA, thf_PA] = changePericenterArg(a,e,omf_PA,om_f,mu);
deltat_CPA_PA_1 = TOF (a, e, theta_PA, thi_PA(1), mu);  % se usi prima intersezione

%fino a pericentro (deltat3)
deltat_BT_PA_1 = TOF (a, e, thf_PA(1), 0, mu);


% Seconda interesezione 
deltat_CPA_PA_2 = TOF (a, e, theta_PA, thi_PA(2), mu);

%fino a pericentro (deltat3)
deltat_BT_PA_2 = TOF (a, e, thf_PA(2), 0, mu);


% Confronto il tempo necessario per compiere le due manovre. 
% NOTA: dalla manovra di cambio di trasferimento in poi la scelta è
% inifluente
if  deltat_BT_PA_1 + deltat_CPA_PA_1 < deltat_BT_PA_2 + deltat_CPA_PA_2 
    y = 1;
    deltat_BT_PA = deltat_BT_PA_1;
    deltat_CPA_PA = deltat_CPA_PA_1;
    fprintf("Per la manovra standard PA il trasferimento più efficiente è nella prima intersezione")
else
    y = 2;
    deltat_BT_PA = deltat_BT_PA_2;
    deltat_CPA_PA = deltat_CPA_PA_2;
    fprintf("Per la manovra standard PA il trasferimento più efficiente è nella seconda intersezione")
end

rr_CP_PA=zeros(3,N_PA);
vv_CP_PA=zeros(3,N_PA);
for j=1:length(th)
[rr_CP_PA(:,j),vv_CP_PA(:,j)]=par2car(a,e,i_f,OM_f,om_f,th(j),mu);
end

th2_PA=(theta_PA+2*pi):dth:(thi_PA(y)+2*pi);
N2_PA=length(th2_PA);
rr_CO1_PA=zeros(3,N2_PA);
vv_CO1_PA=zeros(3,N2_PA);
for j=1:length(th2_PA)
[rr_CO1_PA(:,j),vv_CO1_PA(:,j)]=par2car(a,e,i_f,OM_f,omf_PA,th2_PA(j),mu);
end

[r_CPA_PA,~] = par2car(a, e, i_f, OM_f, omf_PA, thi_PA(y), mu); %manovra nella prima intersezione

figure(Name='manovra cambio di pericentro')
hold on
plot3(rr_CO_PA(1,:),rr_CO_PA(2,:),rr_CO_PA(3,:),'Color',[0.6350, 0.0780, 0.1840])
plot3(rr_CP_PA(1,:),rr_CP_PA(2,:),rr_CP_PA(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(r_COP_PA(1), r_COP_PA(2), r_COP_PA(3), 'o', 'MarkerSize',6,'Color',[0.6350, 0.0780, 0.1840]);
plot3(r_CPA_PA(1), r_CPA_PA(2), r_CPA_PA(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.6350, 0.0780, 0.1840]);
plot3(rr_CO1_PA(1,:),rr_CO1_PA(2,:),rr_CO1_PA(3,:),'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('Orbita cambio di piano','Orbita cambio pericentro','punto di manovra COP','punto di manovra CP')
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)

th3_PA=(thf_PA(y)-2*pi):dth:(2*pi);
N3_PA=length(th3_PA);
rr_CP_PA=zeros(3,N3_PA);
vv_CP_PA=zeros(3,N3_PA);
for j=1:length(th3_PA)
[rr_CP_PA(:,j),vv_CP_PA(:,j)]=par2car(a,e,i_f,OM_f,om_f,th3_PA(j),mu);
end

th3_PA_=(thf_PA(y)):dth:(2*pi);
N3_PA_=length(th3_PA_);
rr_CP_PA_=zeros(3,N3_PA_);
vv_CP_PA_=zeros(3,N3_PA_);
for j=1:length(th3_PA_)
[rr_CP_PA_(:,j),vv_CP_PA_(:,j)]=par2car(a,e,i_f,OM_f,om_f,th3_PA_(j),mu);
end

[r_BT_PA,~] = par2car(a, e, i_f, OM_f, om_f, 0, mu); % inizio trasferimento

figure(Name='fino al pericentro')           
hold on
plot3(rr_CP_PA(1,:),rr_CP_PA(2,:),rr_CP_PA(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(r_CPA_PA(1), r_CPA_PA(2), r_CPA_PA(3), 'o', 'MarkerSize', 6,'Color',[0.9290, 0.6940, 0.1250]);
plot3(r_BT_PA(1), r_BT_PA(2), r_BT_PA(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
plot3(rr_CP_PA_(1,:),rr_CP_PA_(2,:),rr_CP_PA_(3,:),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)
legend('orbita cambio di pericentro','punto di manovra CP','punto di manovra BT1')


%------------------- bitangente (deltat4)

[DeltaV1_BT_PA, DeltaV2_BT_PA,Deltat_PA] = bitangentTransfer (a,e,a_f,e_f,'pa',mu);

rr_BT_PA=zeros(3,N_PA);
vv_BT_PA=zeros(3,N_PA);
for j=1:length(th)
[rr_BT_PA(:,j),vv_BT_PA(:,j)]=par2car(a_f,e_f,i_f,OM_f,om_f,th(j),mu);
end

th5_PA=(pi):dth:(th_f);
N5_PA=length(th5_PA);
rr_BT1_PA=zeros(3,N5_PA);
vv_BT1_PA=zeros(3,N5_PA);
for j=1:length(th5_PA)
[rr_BT1_PA(:,j),vv_BT1_PA(:,j)]=par2car(a_f,e_f,i_f,OM_f,om_f,th5_PA(j),mu);
end

th_T= 0:pi/200:pi;
rp_i=a*(1-e);
ra_f=a_f*(1+e_f);
a_T_PA= (rp_i + ra_f)/2;
e_T_PA= (ra_f - rp_i)/(ra_f +  rp_i);

N_PA=length(th_T);
rr_Tra_PA=zeros(3,N_PA);
vv_Tra_PA=zeros(3,N_PA);
for j=1:length(th)
[rr_Tra_PA(:,j),vv_Tra_PA(:,j)]=par2car(a_T_PA,e_T_PA,i_f,OM_f,om_f,th_T(j),mu);
end

[r_BT2_PA,~] = par2car(a_T_PA, e_T_PA, i_f, OM_f, om_f, pi, mu); % fine trasferimento

figure(Name='trasferimento bitangente PA')          
hold on
plot3(rr_CP_PA(1,:),rr_CP_PA(2,:),rr_CP_PA(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(rr_BT_PA(1,:),rr_BT_PA(2,:),rr_BT_PA(3,:),'Color',[0.4940, 0.1840, 0.5560])
plot3(rr_Tra_PA(1,:),rr_Tra_PA(2,:),rr_Tra_PA(3,:),'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]) 
plot3(r_BT_PA(1), r_BT_PA(2), r_BT_PA(3), 'o', 'MarkerSize', 6,'Color', [0.9290, 0.6940, 0.1250]);
plot3(r_BT2_PA(1), r_BT2_PA(2), r_BT2_PA(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',	[0.4660, 0.6740, 0.1880]);
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('orbita cambio di pericentro','orbita finale','orbita di trasferimento BT','punto di manovra BT1','punto di manovra BT2')
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
exportgraphics(gcf,'confronto.png','Resolution',300)
view(3)

deltat_finale_PA=TOF (a, e, pi, th_f, mu);

[r_BTf_PA,~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu); %punto di arrivo

Delta_man_PA= DeltaV2_BT_PA + DeltaV1_BT_PA +DeltaV_COP_PA + DeltaV_CPA_PA;
delta_t_man_PA = deltat_COP_PA + deltat_BT_PA + deltat_CPA_PA + Deltat_PA + deltat_finale_PA;

figure(Name='fino al punto finale')           
hold on
plot3(rr_BT_PA(1,:),rr_BT_PA(2,:),rr_BT_PA(3,:),'Color',[0.4940, 0.1840, 0.5560])
plot3(r_BT2_PA(1), r_BT2_PA(2), r_BT2_PA(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',	[0.4940, 0.1840, 0.5560]);
plot3(r_BTf_PA(1), r_BTf_PA(2), r_BTf_PA(3), 'o', 'MarkerSize', 6,'Color',[0.4940, 0.1840, 0.5560]);
plot3(rr_BT1_PA(1,:),rr_BT1_PA(2,:),rr_BT1_PA(3,:),'LineWidth',2,'Color',[0.4940, 0.1840, 0.5560])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('orbita finale','punto di manovra BT2','punto di arrivo')
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)


%----------------------- PLOT finale PA 

figure(Name='manovre complete PA')          
hold on
plot3(rr_P_PA(1,:),rr_P_PA(2,:),rr_P_PA(3,:))
plot3(rr_CO_PA(1,:),rr_CO_PA(2,:),rr_CO_PA(3,:))
plot3(rr_CP_PA(1,:),rr_CP_PA(2,:),rr_CP_PA(3,:))
plot3(rr_Tra_PA(1,:),rr_Tra_PA(2,:),rr_Tra_PA(3,:),'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]) 
plot3(rr_BT_PA(1,:),rr_BT_PA(2,:),rr_BT_PA(3,:))

plot3(rr_P1_PA(1,:),rr_P1_PA(2,:),rr_P1_PA(3,:), 'LineWidth',2,'Color',[0, 0.4470, 0.7410])
plot3(rr_CO1_PA(1,:),rr_CO1_PA(2,:),rr_CO1_PA(3,:),'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_CP_PA_(1,:),rr_CP_PA_(2,:),rr_CP_PA_(3,:),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250])
plot3(rr_BT1_PA(1,:),rr_BT1_PA(2,:),rr_BT1_PA(3,:),'LineWidth',2,'Color',[0.4940, 0.1840, 0.5560])

plot3(r_COPi_PA(1), r_COPi_PA(2), r_COPi_PA(3), 'o', 'MarkerSize',6,'Color',	[0, 0.4470, 0.7410]);
plot3(r_COP_PA(1), r_COP_PA(2), r_COP_PA(3), 'o', 'MarkerSize',6,'MarkerFaceColor',	[0, 0.4470, 0.7410]);
plot3(r_CPA_PA(1), r_CPA_PA(2), r_CPA_PA(3), 'o', 'MarkerSize',6,'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
plot3(r_BT_PA(1), r_BT_PA(2), r_BT_PA(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
plot3(r_BT2_PA(1), r_BT2_PA(2), r_BT2_PA(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.4940, 0.1840, 0.5560]);
plot3(r_BTf_PA(1), r_BTf_PA(2), r_BTf_PA(3), 'o', 'MarkerSize', 6,'Color',[0.4940, 0.1840, 0.5560]);
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('Orbita iniziale', 'Orbita cambio di piano', 'Orbita cambio pericentro','Orbita di trasferimento BT PA', 'Orbita finale')
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
camlight
lighting gouraud
exportgraphics(gcf,'orbite complete.png','Resolution',300)
view(45, 30)



%% MANOVRA AP

%-------------------- manovra cambio di piano (deltat1) 
[DeltaV_COP_AP, omf_AP, theta_AP] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);
deltat_COP_AP = TOF (a, e, th_i, theta_AP, mu);

th=th0:dth:thn;
N_AP=length(th);
rr_P_AP=zeros(3,N_PA);
vv_P_AP=zeros(3,N_PA);
for j=1:length(th)
[rr_P_AP(:,j),vv_P_AP(:,j)]=par2car(a,e,i,OM,om,th(j),mu);
end

rr_CO_AP=zeros(3,N_AP);
vv_CO_AP=zeros(3,N_AP);
for j=1:length(th)
[rr_CO_AP(:,j),vv_CO_AP(:,j)]=par2car(a,e,i_f,OM_f,omf_AP,th(j),mu);
end

dth1_AP= (theta_AP+2*pi)/200;
th1_AP=th_i:dth1_AP:(theta_AP+2*pi);
N1_AP=length(th1_AP);
rr_P1_AP=zeros(3,N1_AP);
vv_P1_AP=zeros(3,N1_AP);
for j=1:length(th1_AP)
[rr_P1_AP(:,j),vv_P1_AP(:,j)]=par2car(a,e,i,OM,om,th1_AP(j),mu);
end

[r_COPi_AP,~] = par2car(a, e, i, OM, om, th_i, mu); %punto di inizio
[r_COP_AP,~] = par2car(a, e, i, OM, om, theta_AP, mu); %punto di manovra

figure(Name='manovra cambio di piano') 
hold on
plot3(rr_P_AP(1,:),rr_P_AP(2,:),rr_P_AP(3,:))
plot3(rr_CO_AP(1,:),rr_CO_AP(2,:),rr_CO_AP(3,:))
plot3(r_COPi_AP(1), r_COPi_AP(2), r_COPi_AP(3), 'o', 'MarkerSize',6,'Color',	[0, 0.4470, 0.7410]);
plot3(r_COP_AP(1), r_COP_AP(2), r_COP_AP(3), 'o', 'MarkerSize',6,'MarkerFaceColor',[0, 0.4470, 0.7410]);
plot3(rr_P1_AP(1,:),rr_P1_AP(2,:),rr_P1_AP(3,:), 'LineWidth',2,'Color',[0, 0.4470, 0.7410])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
legend('Orbita iniziale', 'Orbita cambio di piano','punto iniziale','punto di manovra COP')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)


%---------------------- cambio pericentro (deltat2) 

% Prima intersezione 
[DeltaV_CPA_AP, thi_AP, thf_AP] = changePericenterArg(a,e,omf_AP,om_f,mu);
deltat_CPA_AP_1 = TOF (a, e, theta_AP, thi_AP(1), mu);% se usi prima intersezione

%fino a apocentro (deltat3)
deltat_BT_AP_1 = TOF (a, e, thf_AP(1), pi, mu); % da prima intersezione


% Seconda intersezione 
deltat_CPA_AP_2 = TOF (a, e, theta_AP, thi_AP(2), mu);% se usi prima intersezione

%fino a apocentro (deltat3)
deltat_BT_AP_2 = TOF (a, e, thf_AP(2), pi, mu); 


% Confronto il tempo necessario per compiere le due manovre. 
% NOTA: dalla manovra di cambio di trasferimento in poi la scelta è
% inifluente
if  deltat_BT_AP_1 + deltat_CPA_AP_1 < deltat_BT_AP_2 + deltat_CPA_AP_2 
    y = 1;
    deltat_BT_AP = deltat_BT_AP_1;
    deltat_CPA_AP = deltat_CPA_AP_1;
    fprintf("Per la manovra standard AP il trasferimento più efficiente è nella prima intersezione")
else
    y = 2;
    deltat_BT_AP = deltat_BT_AP_2;
    deltat_CPA_AP = deltat_CPA_AP_2;
    fprintf("Per la manovra standard AP il trasferimento più efficiente è nella seconda intersezione")
end

rr_CP_AP=zeros(3,N_AP);
vv_CP_AP=zeros(3,N_AP);
for j=1:length(th)
[rr_CP_AP(:,j),vv_CP_AP(:,j)]=par2car(a,e,i_f,OM_f,om_f,th(j),mu);
end

th2_AP=(theta_AP+2*pi):dth:(thi_AP(y)+2*pi);
N2_AP=length(th2_AP);
rr_CO1_AP=zeros(3,N2_AP);
vv_CO1_AP=zeros(3,N2_AP);
for j=1:length(th2_AP)
[rr_CO1_AP(:,j),vv_CO1_AP(:,j)]=par2car(a,e,i_f,OM_f,omf_AP,th2_AP(j),mu);
end

[r_CPA_AP,~] = par2car(a, e, i_f, OM_f, omf_AP, thi_AP(y), mu); 

th3_AP=(thf_AP(y)-2*pi):dth:(2*pi);
N3_AP=length(th3_AP);
rr_CP1_AP=zeros(3,N3_AP);
vv_CP1_AP=zeros(3,N3_AP);
for j=1:length(th3_AP)
[rr_CP1_AP(:,j),vv_CP1_AP(:,j)]=par2car(a,e,i_f,OM_f,om_f,th3_AP(j),mu);
end

[r_BT1_AP,~] = par2car(a, e, i_f, OM_f, om_f, pi, mu); % inizio trasferimento


figure(Name='manovra cambio di pericentro')
hold on
plot3(rr_CO_AP(1,:),rr_CO_AP(2,:),rr_CO_AP(3,:),'Color',[0.6350, 0.0780, 0.1840])
plot3(rr_CP_AP(1,:),rr_CP_AP(2,:),rr_CP_AP(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(r_COP_AP(1), r_COP_AP(2), r_COP_AP(3), 'o', 'MarkerSize',6,'Color',[0.6350, 0.0780, 0.1840]);
plot3(r_CPA_AP(1), r_CPA_AP(2), r_CPA_AP(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.6350, 0.0780, 0.1840]);
plot3(rr_CO1_AP(1,:),rr_CO1_AP(2,:),rr_CO1_AP(3,:),'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('Orbita cambio di piano','Orbita cambio pericentro','punto di manovra COP','punto di manovra CP')
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)

th3_AP_=(thf_AP(y)-2*pi):dth:(pi);
N3_AP_=length(th3_AP_);
rr_CP1_AP_=zeros(3,N3_AP_);
vv_CP1_AP_=zeros(3,N3_AP_);
for j=1:length(th3_AP_)
[rr_CP1_AP_(:,j),vv_CP1_AP_(:,j)]=par2car(a,e,i_f,OM_f,om_f,th3_AP(j),mu);
end

figure(Name='fino ad apocentro')           
hold on
plot3(rr_CP_AP(1,:),rr_CP_AP(2,:),rr_CP_AP(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(r_CPA_AP(1), r_CPA_AP(2), r_CPA_AP(3), 'o', 'MarkerSize', 6,'Color',[0.9290, 0.6940, 0.1250]);
plot3(r_BT1_AP(1), r_BT1_AP(2), r_BT1_AP(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
plot3(rr_CP1_AP_(1,:),rr_CP1_AP_(2,:),rr_CP1_AP_(3,:),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)
legend('orbita cambio di pericentro','punto di manovra CP','punto di manovra BT1')


%------------------------ bitangente (deltat4) 

[DeltaV1_BT_AP, DeltaV2_BT_AP,Deltat_AP] = bitangentTransfer (a,e,a_f,e_f,'ap',mu);

rr_BT_AP=zeros(3,N_AP);
vv_BT_AP=zeros(3,N_AP);
for j=1:length(th)
[rr_BT_AP(:,j),vv_BT_AP(:,j)]=par2car(a_f,e_f,i_f,OM_f,om_f,th(j),mu);
end

th5_AP=(pi):dth:(th_f);
N5_AP=length(th5_AP);
rr_BT1_AP=zeros(3,N5_AP);
vv_BT1_AP=zeros(3,N5_AP);
for j=1:length(th5_AP)
[rr_BT1_AP(:,j),vv_BT1_AP(:,j)]=par2car(a_f,e_f,i_f,OM_f,om_f,th5_AP(j),mu);
end

th_T= pi:pi/200:2*pi;
ra_t=a*(1+e);
rp_t=a_f*(1-e_f);
a_T_AP= (rp_t + ra_t)/2;
e_T_AP= (ra_t - rp_t)/(ra_t +  rp_t);

N_AP=length(th_T);
rr_Tra_AP=zeros(3,N_AP);
vv_Tra_AP=zeros(3,N_AP);
for j=1:length(th)
[rr_Tra_AP(:,j),vv_Tra_AP(:,j)]=par2car(a_T_AP,e_T_AP,i_f,OM_f,om_f,th_T(j),mu);
end

[r_BT2_AP,~] = par2car(a_T_AP, e_T_AP, i_f, OM_f, om_f, 0, mu); % fine trasferimento


figure(Name='trasferimento bitangente AP')          
hold on
plot3(rr_CP_AP(1,:),rr_CP_AP(2,:),rr_CP_AP(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(rr_BT_AP(1,:),rr_BT_AP(2,:),rr_BT_AP(3,:),'Color',[0.4940, 0.1840, 0.5560])
plot3(rr_Tra_AP(1,:),rr_Tra_AP(2,:),rr_Tra_AP(3,:),'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]) 
plot3(r_BT1_AP(1), r_BT1_AP(2), r_BT1_AP(3), 'o', 'MarkerSize', 6,'Color', [0.9290, 0.6940, 0.1250]);
plot3(r_BT2_AP(1), r_BT2_AP(2), r_BT2_AP(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',	[0.4660, 0.6740, 0.1880]);
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('orbita cambio di pericentro','orbita finale','orbita si trasferimento','punto di manovra BT1','punto di manovra BT2')
axis equal
exportgraphics(gcf,'confronto.png','Resolution',300)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');


deltat_finale_AP=TOF (a, e, 0, th_f, mu);

[r_BTf_PP,~] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu); %punto di arrivo

Delta_man_AP= DeltaV2_BT_AP + DeltaV1_BT_AP +DeltaV_COP_AP + DeltaV_CPA_AP;

delta_t_man_AP = (deltat_COP_AP + deltat_BT_AP + deltat_CPA_AP ...
                 + Deltat_AP + deltat_finale_AP);


%-------------------------- PLOT finale AP 

figure(Name='manovre complete AP')          
hold on
plot3(rr_P_AP(1,:),rr_P_AP(2,:),rr_P_AP(3,:),'Color',	[0, 0.4470, 0.7410])
plot3(rr_CO_AP(1,:),rr_CO_AP(2,:),rr_CO_AP(3,:),'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_CP_AP(1,:),rr_CP_AP(2,:),rr_CP_AP(3,:),'Color',[0.9290, 0.6940, 0.1250])
plot3(rr_Tra_AP(1,:),rr_Tra_AP(2,:),rr_Tra_AP(3,:),'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]) 
plot3(rr_BT_AP(1,:),rr_BT_AP(2,:),rr_BT_AP(3,:),'Color',[0.4940, 0.1840, 0.5560])

plot3(rr_P1_AP(1,:),rr_P1_AP(2,:),rr_P1_AP(3,:), 'LineWidth',2,'Color',[0, 0.4470, 0.7410])
plot3(rr_CO1_AP(1,:),rr_CO1_AP(2,:),rr_CO1_AP(3,:),'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_CP1_AP_(1,:),rr_CP1_AP_(2,:),rr_CP1_AP_(3,:),'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250])

plot3(r_COPi_AP(1), r_COPi_AP(2), r_COPi_AP(3), 'o', 'MarkerSize',6,'Color',	[0, 0.4470, 0.7410]);
plot3(r_COP_AP(1), r_COP_AP(2), r_COP_AP(3), 'o', 'MarkerSize',6,'MarkerFaceColor',	[0, 0.4470, 0.7410]);
plot3(r_CPA_AP(1), r_CPA_AP(2), r_CPA_AP(3), 'o', 'MarkerSize',6,'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
plot3(r_BT1_AP(1), r_BT1_AP(2), r_BT1_AP(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
plot3(r_BT2_AP(1), r_BT2_AP(2), r_BT2_AP(3), 'o', 'MarkerSize', 6,'MarkerFaceColor',[0.4940, 0.1840, 0.5560]);
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
legend('Orbita iniziale', 'Orbita cambio di piano', 'Orbita cambio pericentro','Orbita di trasferimento BT AP', 'Orbita finale')
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
exportgraphics(gcf,'orbite complete.png','Resolution',300)
view(3)

 

%% Confronto strategie standard
% Dati delle strategie
strategie = {'PA', 'AP', ''};
deltaV = [Delta_man_PA,  Delta_man_AP]; % km/s
deltaT = [delta_t_man_PA/3600, delta_t_man_AP/3600]; % h

% Crea figura
figure;
yyaxis left
bar(0.9:1.9 , deltaV, 0.2, 'FaceColor', [0.2 0.6 0.9]);
ylabel('ΔV totale (km/s)')
ylim([0 15])

yyaxis right
bar(1.1:2.1, deltaT, 0.2, 'FaceColor', [0.9 0.4 0.4]);
ylabel('Δt totale (h)')
ylim([0 40])

% Asse X
set(gca, 'xtick', 1:2, 'xticklabel', strategie)
title('Confronto strategie: ΔV e Δt')
exportgraphics(gcf,'confronto.png','Resolution',300)
grid on
 




%% Verifica punto di arrivo dove sta

R_OM_T=[cos(OM_f), sin(OM_f),0;
       -sin(OM_f), cos(OM_f), 0;
       0, 0, 1]; %matrice di rotazione di OM introno a K
R_i_T=[1, 0,0;
       0, cos(i_f), sin(i_f);
       0, -sin(i_f), cos(i_f)]; %matrice di rotazione di i introno a I'
R_om_T=[cos(om_f), sin(om_f),0;
       -sin(om_f), cos(om_f), 0; %matrice di rotazione di i introno a K''
       0, 0, 1];
R_T=(R_om_T*R_i_T*R_OM_T)';

v_per = R_T'*v_f;
r_per = R_T'*r_f;