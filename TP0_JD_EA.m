%% TRABAJO PRÁCTICO TP0 DIMITRI JULIÁN & AUTHIER ELÍAS (2nd Quarter 2024 - UNSAM)

%% Punto 4
close all;
clear;
clc;
%% Item a)
n = 1000000;
X1 = rand(n,1);
X1_max = max(X1);
X1_min = min(X1);
nbins = 1000;
dX1 = (X1_max - X1_min)/nbins;
%% Item b)
h = hist(X1, nbins);
fX1e = h/(dX1*length(X1));
%% Item c)
x1 = X1_min:dX1:X1_max - X1_min;
figure;
plot(x1,fX1e, '-r');
xlabel('X1');
ylabel('f_{X1e}');
title('Gráfico I');
axis([0 1 0.5 1.5]); % Ajustar los límites del eje x y del eje y
%% Item d)
hold on;
fx1 = 1;
plot(x1,fx1*ones(size(x1)), '-b');

%% Item e)
Y1 = X1*pi*3/2 - pi/2;
fy1 = 2/(3*pi);

figure;
plot(Y1, fy1*ones(size(Y1)), '-b');
hold on;
xlabel('Y1');
ylabel('f_{Y1e}');
title('Gráfico II');
axis([-pi/3 pi 0 0.5]); % Ajustar los límites del eje x y del eje y
Y1_e = -pi/2 + (pi + pi/2) * rand(n, 1);
Y1_max = max(Y1_e);
Y1_min = min(Y1_e);
dY1 = (Y1_max - Y1_min)/nbins;
h_y = hist(Y1_e, nbins);
fY1e = h_y/(dY1*length(Y1_e));
y1 = linspace(Y1_min, Y1_max, nbins);
plot(y1,fY1e, '-y');

%% Item f)

figure;
X2 = Y1;
Y2 = sin(X2);
X2_e = Y1_e;
Y2_e = sin(X2_e);
Y2_max = max(Y2_e);
Y2_min = min(Y2_e);
dY2 = (Y2_max - Y2_min)/nbins;
h_y2 = hist(Y2, nbins);
fY2e = h_y2/(dY2*length(Y2_e));
y2a = -0.99:0.01:0;
y2b = 0:0.01:0.99;
fy2a = 2./(pi*3*sqrt(1-y2a.^2));
fy2b = 4./(pi*3*sqrt(1-y2b.^2));
y2 = linspace(Y2_min, Y2_max, nbins);
plot(y2, fY2e,'y');
hold on;
plot(y2a, fy2a, 'r');
hold on;
plot(y2b, fy2b,'r');
axis([-1.5 1.5 0 1]);
xlabel('Y2');
ylabel('f_{Y2}');
title('Gráfico III');
%% Item g)
esperanza_X1= mean(X1)
esperanza_Y1 = mean(Y1)
varianza_X1= var(X1)
varianza_Y1= var(Y1)

%% Punto 6
close all;
clear;
clc;
N = 10000000;
%% Item a)
maxX1 = 1100;
minX1 = 1000;
X1 = minX1 + (maxX1-minX1) .* rand(N,1);

maxX2 = 600;
minX2 = 500;
X2 = minX2 + (maxX2-minX2) .* rand(N,1);

X = [X1,X2];

%% Item b)

nbinsx1 = 75;
nbinsx2 = 75;
X1_max = max(X1);
X1_min = min(X1);
dX1 = (X1_max - X1_min)/nbinsx1;
X2_max = max(X2);
X2_min = min(X2);
dX2 = (X2_max - X2_min)/nbinsx2;

fX1X2e = hist3(X, [nbinsx1 nbinsx2])/(dX1*dX2*N);
fX1X2 = 1/10000 * ones([nbinsx1, nbinsx2]);

%% Item c)
x1 = X1_min:dX1:X1_max - dX1;
x2 = X2_min:dX2:X2_max - dX2;
[x1,x2] = meshgrid(x1, x2);
mesh(x1, x2, fX1X2e);
hold on;
mesh(x1, x2, fX1X2,'EdgeColor', 'r');

%% Item d)
M = [[3,2];[2,2]];
Y = M*X';
Y1 = Y(1,:);
Y2 = Y(2,:);

nbinsy1 = 75;
nbinsy2 = 75;
Y1_max = max(Y1);
Y1_min = min(Y1);
dY1 = (Y1_max - Y1_min)/nbinsy1;
Y2_max = max(Y2);
Y2_min = min(Y2);
dY2 = (Y2_max - Y2_min)/nbinsy2;

figure;
fY1Y2e = hist3(Y', [nbinsy1 nbinsy2])/(dY1*dY2*N);

y1 = Y1_min:dY1:Y1_max - dY1;
y2 = Y2_min:dY2:Y2_max - dY2;
[y1,y2] = meshgrid(y1, y2);
mesh(y1, y2, fY1Y2e);


%% Item e)
figure;
contour(y1,y2,fY1Y2e);
axis([3900 4600 2900 3500]);

%% Item f)
%P(4300<Y1<4500)
P1 = length(find( (Y1>4300) & (Y1<4500)))/N
%P(3000<Y2<3200)
P2 = length(find( (Y2>3000) & (Y2<3200)))/N

%% Item g)
fX1 = 1/100;
fX2 = 1/100;
Y1 = 3.*X1 + 2.*X2;
Y2 = 2.*X1 + 2.*X2;

%Esperanzas
esp_X1=mean(X1)
esp_X2=mean(X2)
esp_X1X2=mean(X1.*X2)
esp_X1X1=mean(X1.^2)
esp_X2X2=mean(X2.^2)
esp_Y1=mean(Y1)
esp_Y2=mean(Y2)
esp_Y1Y2=mean(Y1.*Y2)
esp_Y1Y1=mean(Y1.^2)
esp_Y2Y2=mean(Y2.^2)


%% Item h)
var_X1=var(X1)
var_X2=var(X2)
cov_X1X2_matriz=cov(X1,X2)
cov_X1X2 = cov_X1X2_matriz(1,2)
coef_corrX1X2= cov_X1X2/((var_X1*var_X2)^0.5)
var_Y1=var(Y1)
var_Y2=var(Y2)
cov_Y1Y2_matriz=cov(Y1,Y2)
cov_Y1Y2=cov_Y1Y2_matriz(1,2)
coef_corrY1Y2= cov_Y1Y2/((var_Y1*var_Y2)^0.5)

