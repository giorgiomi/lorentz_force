%PAF script
clear; close all; clc;
addpath scripts/;
g = 9.806;

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

fit_data = regressione_lineare(F,V,dV)

xx = (0:0.1e-03:25e-03)*g;
plot(xx,fit_data.m*xx+fit_data.b);

fatt_conv = fit_data.m; %in V/N
dfatt_conv = fit_data.dm;
fprintf("Fattore di conversione: %.4f +/- %.5f V/N\n",fatt_conv,dfatt_conv);

%data from DMM
%{
dataDMM = readmatrix("defbuffer1.csv");
figure();
plot(dataDMM(:,1),dataDMM(:,2),'.');
%}
