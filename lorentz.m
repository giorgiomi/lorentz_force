%% Lorentz force measure
clear; close all; clc;
addpath scripts/;
g = 9.806;

%% Calibrazione sensore di forza
RSG_verde_giallo = 0.99855e03;
R_rosso_nero = 1.05177e03;
V_BIAS = 12;
i_rosso_nero = V_BIAS/R_rosso_nero;
i_rosso_nero_meas = 11e-03; %per ora le cose funzionano.
V_offset = 432e-06; %432 uV

data = readmatrix("data.xlsx");
m = data(1:end,1); F = m*g;
V = data(1:end,2);
dV = data(:,3);

figure();
axes();
hold on
grid on
errorbar(F,V,dV,'.','MarkerSize',12);
xlabel("force [N]");
ylabel("voltage [V]")

fit_data = regressione_lineare(F,V,dV);

xx = (0:0.1e-03:25e-03)*g;
plot(xx,fit_data.m*xx+fit_data.b);

fatt_conv = fit_data.m; %in V/N
dfatt_conv = fit_data.dm;
fprintf("Fattore di conversione: %.4f +/- %.5f V/N\n",fatt_conv,dfatt_conv);


%% Accoppiamento in DC
rr = 4.65; %resistenza bobina ricevitrice in ohm
rs = 4.64; %resistenza bobina sorgente in ohm

dataV = readmatrix("voltage_sensor.csv");
dataI = readmatrix("current_sorg.csv");
figure();
hold on
grid on
a = 2600;
b = 5500;
t = dataV(:,1); t = t(a:b);
V = dataV(:,2); V = V(a:b);
tI = dataI(:,1); tI = tI(2.34e4:end);
I = dataI(:,2); I = I(2.34e4:end);
plot(t,V,'.-');

figure();
hold on;
grid on;
plot(tI,I,'.-');


%% Studio dalla dipendenza da VDC
d = 1.115e-02; %in cm
VDC = [6.113 5.114 4.112 3.110 2.106 2.609 3.612];










