clear
clc
close all

% Dati asteroide 
a_A=1.430330;              % UA
a_A=a_A*1.496e8;           % Km
e_A= 0.256103;
i_A= 8.71;                 % gradi
i_A= i_A*2*pi/360;         % rad
OM_A= 246.29;              % gradi
OM_A= OM_A *2*pi/360;      % rad
om_A= 338.44;              % gradi
om_A= om_A*2*pi/360;       % rad
mu_S= 132712440018;

% Dati Terra
a_T= 1.4946e8;
e_T= 0.016;
i_T= 9.1920e-5;
OM_T= 2.7847;
om_T= 5.2643;
mu_T=398600;

% Plot orbite della terra e dell'asteroide
th0=0;
thn=2*pi;
dth= 2*pi/200;
th=th0:dth:thn;

N=length(th);
rr_T=zeros(3,N);
vv_T=zeros(3,N);
for j=1:length(th)
[rr_T(:,j),vv_T(:,j)]=par2car(a_T,e_T,i_T,OM_T,om_T,th(j),mu_S);
end

rr_A=zeros(3,N);
vv_A=zeros(3,N);

for j=1:length(th)
[rr_A(:,j),vv_A(:,j)]=par2car(a_A,e_A,i_A,OM_A,om_A,th(j),mu_S);
end

%scelgo 2 parametri da fissare
%inclinazione fissata come la terra da cui derivo theta finale (una tra le 
% 2 intersezioni dell'orbita finale col piano della terra) e ilo punto di
% partenza nel pericentro dell'orbita iniziale
i_tra=i_T;
th_T1=0;

%per i punti finali trovo l'intersezione tra i piani orbitali (linea dei nodi)
R_OM_T=[cos(OM_T), sin(OM_T),0;
       -sin(OM_T), cos(OM_T), 0;
       0, 0, 1]; %matrice di rotazione di OM introno a K
R_i_T=[1, 0,0;
       0, cos(i_T), sin(i_T);
       0, -sin(i_T), cos(i_T)]; %matrice di rotazione di i introno a I'
R_om_T=[cos(om_T), sin(om_T),0;
       -sin(om_T), cos(om_T), 0; %matrice di rotazione di i introno a K''
       0, 0, 1];
R_T=(R_om_T*R_i_T*R_OM_T)';

R_OM_A=[cos(OM_A), sin(OM_A),0;
       -sin(OM_A), cos(OM_A), 0;
       0, 0, 1]; 
R_i_A=[1, 0,0;
       0, cos(i_A), sin(i_A);
       0, -sin(i_A), cos(i_A)];
R_om_A=[cos(om_A), sin(om_A),0;
       -sin(om_A), cos(om_A), 0; 
       0, 0, 1];
R_A=(R_om_A*R_i_A*R_OM_A)';

n_T=R_T(:,3); %versore normale piano della terra
n_A=R_A(:,3); %versore nromale piano asteroide
p_T=R_T(:,1);
p_A=R_A(:,1); %direzione pericentro orbita asteroide

N = (cross(n_T, n_A))/norm(cross(n_T, n_A)); %linea dei nodi

th_asc = acos(dot(N, p_A)); %nodo ascendente
th_desc = th_asc+pi; %nodo discendente

[rr_A2_,vv_A2_]=par2car(a_A,e_A,i_A,OM_A,om_A,th_asc,mu_S);
[rr_A2__,vv_A2__]=par2car(a_A,e_A,i_A,OM_A,om_A,th_desc,mu_S);

figure
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:),'Color',[0, 0.4470, 0.7410])
hold on
plot3(rr_A(1,:),rr_A(2,:),rr_A(3,:),'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_A2_(1), rr_A2_(2), rr_A2_(3), 'o', 'MarkerSize',6,'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
plot3(rr_A2__(1), rr_A2__(2), rr_A2__(3), 'o', 'MarkerSize',6,'MarkerFaceColor', [0.4940, 0.1840, 0.5560]);

opts.Position = [0,0,0];
opts.Units = 'ft';
planet3D('Sun',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
view(3)
exportgraphics(gcf,'orbita iniziale.png','Resolution',300)
legend('Orbita Terra','Orbita asteroide','nodo ascendente','nodo discendente')


%% nodo ascendente
th_A2=th_asc;

%definisco i versori
i=[1;0;0];
j=[0;1;0];
k=[0;0;1];

%RAAN
OM_tra=OM_T; %(stessa inclinazione --> stessa RAAN)

%vettore di stato nel punto di partenza(Terra) 
[rr_T1,vv_T1]=par2car(a_T,e_T,i_T,OM_T,om_T,th_T1,mu_S);
    
%vettore di stato nel punto di arrivo(asteroide) 
[rr_A2,vv_A2]=par2car(a_A,e_A,i_A,OM_A,om_A,th_A2,mu_S);
    
% Matrice di trasferimento
  R_OM_tra=[cos(OM_tra), sin(OM_tra),0;
            -sin(OM_tra), cos(OM_tra), 0;
             0, 0, 1]; 
  R_i_tra=[1, 0,0;
           0, cos(i_tra), sin(i_tra);
           0, -sin(i_tra), cos(i_tra)];

Delta_V_min_asc = inf;

for om_tra=0:pi/200:2*pi %scelgo la variabile libera su cui far ciclare (per minimizzare deltav)
    R_om_tra=[cos(om_tra), sin(om_tra),0;
             -sin(om_tra), cos(om_tra), 0; 
              0, 0, 1];
    R_tra=(R_om_tra*R_i_tra*R_OM_tra)';
    %R_tra: trasforma da sistema locale (orbita) a inerziale
    %R_tra': trasforma da inerziale a sistema orbitale (local frame)
    rr_tra1=R_tra'*rr_T1; %posizione di partenza nel sdr dell'orbita di trasferimento
    rr_tra2=R_tra'*rr_A2; %posizione di arrivo nel sdr dell'orbita di trasferimento
   
    theta1 = atan2(rr_tra1(2), rr_tra1(1));
    if theta1 < 0
    theta1 = theta1 + 2*pi;
    end
    theta2 = atan2(rr_tra2(2), rr_tra2(1));
    if theta2 < 0
    theta2 = theta2 + 2*pi;
    end

    r1 = norm(rr_tra1);
    r2 = norm(rr_tra2);

    e_tra = (r2 - r1) / (r1 * cos(theta1) - r2 * cos(theta2));

    if e_tra < 0 || e_tra >= 1 || ~isreal(e_tra)
    continue  % salta orbite aperte, paraboliche o non fisiche
    end

    p_tra = r1 * (1 + e_tra * cos(theta1));
    a_tra = p_tra / (1 - e_tra^2);

    [~, v_tra1] = par2car (a_tra, e_tra, i_tra, OM_tra, om_tra, theta1, mu_S);
    [~, v_tra2] = par2car (a_tra, e_tra, i_tra, OM_tra, om_tra, theta2, mu_S);

    delta_V_1 = norm(v_tra1-vv_T1);
    delta_V_2 = norm(vv_A2-v_tra2);
    
   if delta_V_2 + delta_V_1 < Delta_V_min_asc
        Delta_V_min_asc = delta_V_2 + delta_V_1;
        best_om_asc = om_tra;
        best_e_asc = e_tra;
        best_a_asc = a_tra;
        best_theta1_asc = theta1;
        best_theta2_asc = theta2;
    end
end

[best_r1_asc, best_v1_asc] = par2car(best_a_asc, best_e_asc, i_tra, OM_tra, best_om_asc, best_theta1_asc, mu_S);
[best_r2_asc, best_v2_asc] = par2car(best_a_asc, best_e_asc, i_tra, OM_tra, best_om_asc, best_theta2_asc, mu_S);
[rr_T1_]=par2car(a_T,e_T,i_T,OM_T,om_T,th_T1,mu_S);
[rr_A2_,vv_A2_]=par2car(a_A,e_A,i_A,OM_A,om_A,th_asc,mu_S);

tof_asc=TOF(best_a_asc,best_e_asc,best_theta1_asc,best_theta2_asc,mu_S);

th_T=linspace(best_theta2_asc+2*pi,best_theta1_asc);
N=length(th_T);
rr_Tra=zeros(3,N);
vv_Tra=zeros(3,N);
for j=1:N
    [rr_tra_opt(:,j),~]=par2car(best_a_asc,best_e_asc,i_tra,OM_tra,best_om_asc, th_T(j), mu_S);
end

figure        
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:),'Color',[0, 0.4470, 0.7410])
hold on
plot3(rr_A(1,:),rr_A(2,:),rr_A(3,:),'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_T1_(1), rr_T1_(2), rr_T1_(3), 'o', 'MarkerSize',6,'MarkerFaceColor',	[0, 0.4470, 0.7410]);
plot3(rr_A2_(1), rr_A2_(2), rr_A2_(3), 'o', 'MarkerSize',6,'Color',[0.8500, 0.3250, 0.0980]);
plot3(rr_tra_opt(1,:), rr_tra_opt(2,:), rr_tra_opt(3,:), '--','LineWidth',1.5,'Color',[0.4660, 0.6740, 0.1880])

opts.Position = [0,0,0];
opts.Units = 'ft';
planet3D('Sun',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('Orbita Terra','Orbita Asteroide', 'Punto di partenza','Punto di arrivo', 'Orbita di trasferimento')
exportgraphics(gcf,'orbita iniziale.png','Resolution',300)

view(3)

%% nodo discendente
th_A2=th_desc;

%definisco i versori
i=[1;0;0];
j=[0;1;0];
k=[0;0;1];

%RAAN
OM_tra=OM_T; %(stessa inclinazione --> stessa RAAN)

%vettore di stato nel punto di partenza(Terra) 
[rr_T1,vv_T1]=par2car(a_T,e_T,i_T,OM_T,om_T,th_T1,mu_S);
    
%vettore di stato nel punto di arrivo(asteroide) 
[rr_A2,vv_A2]=par2car(a_A,e_A,i_A,OM_A,om_A,th_A2,mu_S);
    
% Matrice di trasferimento
  R_OM_tra=[cos(OM_tra), sin(OM_tra),0;
            -sin(OM_tra), cos(OM_tra), 0;
             0, 0, 1]; 
  R_i_tra=[1, 0,0;
           0, cos(i_tra), sin(i_tra);
           0, -sin(i_tra), cos(i_tra)];

Delta_V_min_desc = inf;

for om_tra=0:pi/200:2*pi

    R_om_tra=[cos(om_tra), sin(om_tra),0;
             -sin(om_tra), cos(om_tra), 0; 
              0, 0, 1];

    R_tra=(R_om_tra*R_i_tra*R_OM_tra)';
    %R_tra: trasforma da sistema locale (orbita) a inerziale
    %R_tra': trasforma da inerziale a sistema orbitale (local frame)
    rr_tra1=R_tra'*rr_T1; %posizione di partenza nel sdr dell'orbita di trasferimento
    rr_tra2=R_tra'*rr_A2; %posizione di arrivo nel sdr dell'orbita di trasferimento

    theta1 = atan2(rr_tra1(2), rr_tra1(1));
    if theta1 < 0
    theta1 = theta1 + 2*pi;
    end
    theta2 = atan2(rr_tra2(2), rr_tra2(1));
    if theta2 < 0
    theta2 = theta2 + 2*pi;
    end

    r1 = norm(rr_tra1);
    r2 = norm(rr_tra2);

    e_tra = (r2 - r1) / (r1 * cos(theta1) - r2 * cos(theta2));

    if e_tra < 0 || e_tra >= 1 || ~isreal(e_tra)
    continue  % salta orbite aperte, paraboliche o non fisiche
    end

    p_tra = r1 * (1 + e_tra * cos(theta1));
    a_tra = p_tra / (1 - e_tra^2);

    [~, v_tra1] = par2car (a_tra, e_tra, i_tra, OM_tra, om_tra, theta1, mu_S);
    [~, v_tra2] = par2car (a_tra, e_tra, i_tra, OM_tra, om_tra, theta2, mu_S);

    delta_V_1 = norm(v_tra1-vv_T1);
    delta_V_2 = norm(vv_A2-v_tra2);
    
   if delta_V_2 + delta_V_1 < Delta_V_min_desc
        Delta_V_min_desc = delta_V_2 + delta_V_1;
        best_om_desc = om_tra;
        best_e_desc = e_tra;
        best_a_desc = a_tra;
        best_theta1_desc = theta1;
        best_theta2_desc = theta2;
    end
end

[best_r1_desc, best_v1_desc] = par2car(best_a_desc, best_e_desc, i_tra, OM_tra, best_om_desc, best_theta1_desc, mu_S);
[best_r2_desc, best_v2_desc] = par2car(best_a_desc, best_e_desc, i_tra, OM_tra, best_om_desc, best_theta2_desc, mu_S);
[rr_T1_]=par2car(a_T,e_T,i_T,OM_T,om_T,th_T1,mu_S);
[rr_A2_,vv_A2_]=par2car(a_A,e_A,i_A,OM_A,om_A,th_asc,mu_S);

tof_desc=TOF(best_a_desc,best_e_desc,best_theta1_desc,best_theta2_desc,mu_S);

th_T=linspace(best_theta2_desc+2*pi,best_theta1_desc);
N=length(th_T);
rr_Tra=zeros(3,N);
vv_Tra=zeros(3,N);
for j=1:N
    [rr_tra_opt(:,j),~]=par2car(best_a_desc,best_e_desc,i_tra,OM_tra,best_om_desc, th_T(j), mu_S);
end

fig = figure;        
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:),'Color',[0, 0.4470, 0.7410])
hold on
plot3(rr_A(1,:),rr_A(2,:),rr_A(3,:),'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_T1_(1), rr_T1_(2), rr_T1_(3), 'o', 'MarkerSize',6,'MarkerFaceColor',	[0, 0.4470, 0.7410]);
plot3(rr_A2__(1), rr_A2__(2), rr_A2__(3), 'o', 'MarkerSize',6,'Color',[0.8500, 0.3250, 0.0980]);
plot3(rr_tra_opt(1,:), rr_tra_opt(2,:), rr_tra_opt(3,:), '--','LineWidth',1.5,'Color',[0.4660, 0.6740, 0.1880])

opts.Position = [0,0,0];
opts.Units = 'ft';
planet3D('Sun',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('Orbita Terra','Orbita Asteroide', 'Punto di partenza','Punto di arrivo', 'Orbita di trasferimento')
exportgraphics(gcf,'orbita iniziale.png','Resolution',300)

view(3)

%% verifica manovra migliore

if Delta_V_min_asc < Delta_V_min_desc
    disp('Scelto nodo ascendente');
    best_r1 = best_r1_asc;
    best_r2 = best_r2_asc;
    best_v1 = best_v1_asc;
    best_v2 = best_v2_asc;
    Delta_V_min = Delta_V_min_asc;
else
    disp('Scelto nodo discendente');
    best_r1 = best_r1_desc;
    best_r2 = best_r2_desc;
    best_v1 = best_v1_desc;
    best_v2 = best_v2_desc;
    Delta_V_min = Delta_V_min_desc;
end








