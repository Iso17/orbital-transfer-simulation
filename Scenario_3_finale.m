clear
clc
close all

%orbita di parcheggio satellite @Terra
r_T=[-1524.866
     9996.648	
     4556.442]; 
v_T=[-6.421	
     -2.655
     3.6760];
mu_T=398600;
[a_T,e_T,i_T,OM_T,om_T,th_T] = car2par(r_T, v_T, mu_T);
M_T=5.974e+24;

%orbita di trasferimento del satellite @Sole
mu_S=132712440018;
a_S=1.5524e+08;
e_S=0.0527;
i_S=9.1920e-05;
OM_S=2.7847;
om_S=5.3250;
theta1=6.2225;
theta2=2.4727;
M_S=1.989e+30;

%orbita asteroide @Sole (non rilevante per le patched chonics)
a_A=1.430330;              % UA
a_A=a_A*1.496e8;           % Km
e_A= 0.256103;
i_A= 8.71;                 % gradi
i_A= i_A*2*pi/360;         % rad
OM_A= 246.29;              % gradi
OM_A= OM_A *2*pi/360;      % rad
om_A= 338.44;              % gradi
om_A= om_A*2*pi/360;       % rad
M_A=1.2741e+13; 
mu_A=M_A*6.67e-11;
r_A=1.15;

%velocità utili
v_tra1=[-30.0681
  -11.1830
     0.0019]; %velocità del satellite nell'orbita di trasferimento nel punto di partenza @sole
v_tra2=[23.9885
  -17.1041
    0.0007]; %velocità del satellite nell'orbita di trasferimento nel punto di arrivo @sole
v_1i=[-29.7051
  -5.8667
    0.0015]; %velocità terra nel punto di partenza @sole
v_2f=[27.9017
  -14.8849
    4.8307]; %velocità asteroide nel punto di arrivo @sole

%plot delle tre orbite da collegare
th0=0;
thn=2*pi;
dth= 2*pi/200;
th=th0:dth:thn;
N=length(th);

rr_S=zeros(3,N);
vv_S=zeros(3,N);
for j=1:length(th)
[rr_S(:,j),vv_S(:,j)]=par2car(a_S,e_S,i_S,OM_S,om_S,th(j),mu_S);
end

rr_A=zeros(3,N);
vv_A=zeros(3,N);
for j=1:length(th)
[rr_A(:,j),vv_A(:,j)]=par2car(a_A,e_A,i_A,OM_A,om_A,th(j),mu_S);
end

figure(Name='orbite eliocentriche')  
%l'orbita di parcheggio non può apparire perchè troppo piccola rispettlo le
%orbite eliocentriche, serve un grafico a parte
hold on
plot3(rr_S(1,:),rr_S(2,:),rr_S(3,:),'Color',[0.8500, 0.3250, 0.0980])
plot3(rr_A(1,:),rr_A(2,:),rr_A(3,:),'Color',[0.9290, 0.6940, 0.1250])

opts.Position = [0,0,0];
opts.Units = 'ft';
planet3D('Sun',opts)

opts.Position = [149597870,0,0 ];
opts.Units = 'm';
planet3D('Earth',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('orbita eliocentrica','orbita asteroide');
  

rr_T=zeros(3,N);
vv_T=zeros(3,N);
for j=1:length(th)
[rr_T(:,j),vv_T(:,j)]=par2car(a_T,e_T,i_T,OM_T,om_T,th(j),mu_T);
end

figure(Name='orbita geocentrica')  
hold on
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:),'Color',[0, 0.4470, 0.7410])
opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

%--------------------------------------------------------------------------
%prima iperbole: sfuggire dall'attrazione della terra
v_inf1=norm(v_tra1-v_1i);%velocità di eccesso iperbolico-perbole di fuga

r_op = a_T * (1 - e_T); %raggio di pericetro=raggio di pericentro orbita di parcehggio
v_op = sqrt(mu_T * (2/r_op - 1/a_T)); %velocità dell'orbita di parcheggio (non è circolare)

a_SOI_T = 151590000; %distanza terra-sole
r_SOI_T = a_SOI_T * (M_T / M_S)^(2/5);

a_1=-mu_T/((v_inf1)^2);
e_1=1+r_op*v_inf1^2/mu_T;
v_f1=sqrt(2*mu_T/r_op); %velocità di fuga
v_ph1=sqrt(v_inf1^2+v_f1^2); %velocità al pericentro dell'iperbole

delta_n1=-a_1*sqrt(e_1^2-1);%parametro d'impatto
delta_1=2 * asin(1/e_1);%angolo di deflessione

%deltav per inserire il satellite nella traiettoria di fuga
delta_v1=v_ph1-v_op;

%plot iperbole di fuga
theta_h = linspace(0, pi/1.7, 500);
r_1 = zeros(3, length(theta_h));
v_1 = zeros(3, length(theta_h));
for j = 1:length(theta_h)
    [r_1(:,j), v_1(:,j)] = par2car(a_1, e_1, i_T, OM_T, om_T, theta_h(j), mu_T);
end
theta_h_ = linspace(-pi/2, 0, 500);
r_1_ = zeros(3, length(theta_h_));
v_1_ = zeros(3, length(theta_h_));
for j = 1:length(theta_h)
    [r_1_(:,j), v_1_(:,j)] = par2car(a_1, e_1, i_T, OM_T, om_T, theta_h_(j), mu_T);
end

figure
hold on
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:),'Color',[0, 0.4470, 0.7410]);
plot3(r_1(1,:),r_1(2,:),r_1(3,:),'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
plot3(r_1_(1,:),r_1_(2,:),r_1_(3,:),'--','Color',[0.8500, 0.3250, 0.0980]);
plot3(rr_T(1),rr_T(2),rr_T(3),'o', 'MarkerSize',6,'MarkerFaceColor',[0.8500, 0.3250, 0.0980])

opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
legend('Orbita di parcheggio','Iperbole di fuga')
exportgraphics(gcf,'iperbole di fuga.png','Resolution',300)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');


%--------------------------------------------------------------------------
%seconda iperbole: entrare nell'orbita di cattura dell'asteroide
v_inf2=norm(v_2f-v_tra2);%velocità di eccesso iperbolico-perbole di arrivo

a_SOI_A = a_A; %distanza media asteroide-sole
r_SOI_A = a_SOI_A * (M_A / M_S)^(2/5); %sphere of influence

a_2 = -mu_A / v_inf2^2;
 
r_ph2_max = r_SOI_A;  
r_ph2=linspace(r_A, r_ph2_max, 5);

e_2 = zeros(1, 5);
delta_v2 = zeros(1, 5);
delta_2 = zeros(1, 5);
delta_n2 = zeros(1, 5);  
r_oc = zeros(1,5);
v_ph2 = zeros(1,5);
v_oc=zeros(1,5);
v_f2=zeros(1,5);

for i = 1:5
    r_oc(i) = r_ph2(i); 
    e_2(i) = 1 + r_oc(i) * v_inf2^2 / mu_A;
    v_f2(i) = sqrt(2 * mu_A / r_oc(i)); %velocità di fuga
    v_oc(i) = sqrt(mu_A/r_oc(i) ); %velocità orbite circolari 
    v_ph2(i) = sqrt(v_inf2^2 + v_f2(i)^2); %velocità al pericentro delle iperboli
    delta_n2(i) = -a_2 * sqrt(e_2(i)^2 - 1);
    delta_2(i) = 2 * asin(1/e_2(i));
    delta_v2(i) = norm(v_oc(i)-v_ph2(i));
end

figure
plot(r_oc, delta_v2, '*')
xlabel('Raggio di pericentro r_{oc} [km]')
ylabel('\Delta v_2 [km/s]')
title('Delta V per inserimento in orbita iperbolica')
grid on

figure
plot(r_oc, delta_n2, '*')
xlabel('Raggio di pericentro r_{oc} [km]')
ylabel('\delta_{n2} [km]')
title('Parametro d''impatto \delta_{n2}')
grid on

%plot iperbole di cattura (per varie orbite circolari attorno
%all'asteroide)

figure
hold on

for i=1:5
    theta_A = linspace(0, 2*pi, 300);
    circle_A = [cos(theta_A); sin(theta_A); zeros(1,length(theta_A))] * r_oc(i);

    plot3(circle_A(1,:), circle_A(2,:), circle_A(3,:), '--', 'HandleVisibility', 'off');

    theta_h2 = linspace(0, pi/2, 300);
    r_h2 = zeros(3, length(theta_h2));
    for j = 1:length(theta_h2)
        [r_h2(:,j), ~] = par2car(a_2, e_2(i), i_A, OM_A, om_A, theta_h2(j), mu_A);
    end

    plot3(r_h2(1,:), r_h2(2,:), r_h2(3,:),'DisplayName', ['iperbole ' num2str(i)]);
end
[X, Y, Z] = sphere(50); 
X = X * r_A;
Y = Y * r_A;
Z = Z * r_A;
surf(X, Y, Z, 'HandleVisibility', 'off');
legend('iperbole 1','iperbole 2','iperbole 3','iperbole 4', 'iperbole 5')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

grid on
axis equal
exportgraphics(gcf,'iperbole di cattura.png','Resolution',300)



%% scenario 3 bis
%iperbole di fuga su un piano diverso rispetto l'orbita di parcheggio
% rotazione del sistema di riferimento da planetario (terra) a ecclittico 
% del sole (ecclittico)

eps = deg2rad(23.45); % inclinazione orbita terra 

T_eci = [1 0 0
         0 cos(eps) sin(eps)
         0 -sin(eps) cos(eps)]; % matrice di rotazione intorno all'asse x
v_inf_ecl = v_tra1 - v_1i; %eccesso iperbolico in ECLIP
v_inf_Eci = T_eci' * v_inf_ecl; %eccesso iperbolico in ECI

r_inf = v_inf_Eci/norm(v_inf_Eci);% direzione del asintoto uscente dell'iperbole 
r_h = norm(r_T); %intersezione tra iperbole e orbita di parcheggio (NON è il pericentro 
% dell'iperbole)
a_H = -mu_T/norm(v_inf_Eci)^2;
alpha = acos(dot(r_T, r_inf)/(norm(r_T)*norm(r_inf)));

% sistema non lineare
% x = [e_H, theta_H, theta_inf]
fun = @(x) [ 
r_h - a_H*(1-x(1).^2)/(1 + x(1)*cos(x(2)));
cos(x(3)) + 1/x(1);
x(3) - alpha - x(2);
];

x_0 = [1.5, 0.1, 0.5];% Guess iniziale 

% Risoluzione del sistema non lineare
options = optimoptions('fsolve','Display','off');
x = fsolve(fun, x_0, options);

e_H = x(1);
theta_H = x(2);
theta_inf = x(3);
if x(3)<pi
    fprintf("soluzione accettabile")
else 
    fprintf("soluzione NON accettabile")
end

r_ph = a_H * (1 - x(1));

if r_ph > 6378 %raggio terra in km 
        fprintf("\n\nsoluzione accettabile")
else 
    fprintf("\n\nsoluzione NON accettabile")
end

fun_2 = @(e_H_vers) [
    dot(r_T, e_H_vers) - cos(theta_H);
    dot(r_inf, e_H_vers) - cos(theta_inf);
    norm(e_H_vers) - 1;
    ];

e_H_0 = [0,0,1];

% Risoluzione del sistema non lineare
options = optimoptions('fsolve','Display','off');
e_H_vers = fsolve(fun_2, e_H_0, options);

%piano su cui giace l'iperbole
h_H = cross(r_T, v_inf_Eci);
h_H = h_H / norm(h_H);  % versore

k = [0, 0, 1]; 
i_H = acos(dot(h_H, k));

N_H = cross(k, h_H); 
N_H = N_H / norm(N_H);
 
i = [1, 0, 0];
if N_H(2) >= 0
    OM_H = acos(dot(N_H, i));
else
    OM_H = 2*pi - acos(dot(N_H, i));
end

cos_beta=dot(N_H,r_T)/(norm(N_H)*norm(r_T));
beta=acos(cos_beta);
om_H=beta-theta_H;

[~,v_h_T] = par2car(a_H, e_H, i_H, OM_H, om_H, theta_H, mu_T); %velocità iperbole nell'intersezione con l'orbita di parcheggio
delta_v_bis = norm(v_h_T - v_T);

%plot iperbole di fuga 
theta_hh1 = linspace(-1.7, theta_H, 1000);
r_h1 = zeros(3, length(theta_hh1));
v_h1 = zeros(3, length(theta_hh1));
for j1 = 1:length(theta_hh1)
    [r_h1(:,j1), v_h1(:,j1)] = par2car(a_H, x(1),i_H ,OM_H ,om_H , theta_hh1(j1), mu_T);
end

theta_hh2 = linspace(theta_H, 1.9, 1000);
r_h2 = zeros(3, length(theta_hh2));
v_h2 = zeros(3, length(theta_hh2));
for j2 = 1:length(theta_hh2)
    [r_h2(:,j2), v_h2(:,j2)] = par2car(a_H, x(1),i_H ,OM_H ,om_H , theta_hh2(j2), mu_T);
end

figure
hold on
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:),'Color',[0, 0.4470, 0.7410]);
plot3(r_h1(1,:),r_h1(2,:),r_h1(3,:),'--','Color',[0.8500, 0.3250, 0.0980]);
plot3(r_h2(1,:),r_h2(2,:),r_h2(3,:),'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
plot3(r_T(1), r_T(2), r_T(3), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);  % Punto di intersezione

opts.Position = [0,0,0];
opts.Units = 'km';
planet3D('Earth',opts)

grid on
axis equal
legend('Orbita di parcheggio','Iperbole di fuga')
exportgraphics(gcf,'iperbole di fuga.png','Resolution',300)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

