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
t = dataV(:,1); t = t(a:b); t = t/2;
V = dataV(:,2); V = V(a:b);
tI = dataI(:,1); tI = tI(2.34e4:end);
I = dataI(:,2); I = I(2.34e4:end);
plot(t,V,'.-');
xlabel("t"); ylabel("V");

%{
figure();
hold on;
grid on;
plot(tI,I,'.-');
%}


%% Fit sinusoidale
close all; clc;
figtimeseries = figure;
plot(t, V); % Assuming it's a temporal series (t, V)
hold on; grid on;
vlim = [0 15e-3];
ylim(vlim);
title("Figura da fittare");

% Define intervals to analyze at different distances
tran = [10000 15000; 15200 20000];

figure(figtimeseries);
for jj = 1:size(tran, 1)
    plot(tran(jj, 1) * [1 1], vlim, 'g--');
    plot(tran(jj, 2) * [1 1], vlim, 'm--');
end

tOFF = [1300; 1300]; % Offset t0, manually set
fmod = 0.1; % Modulation frequency
Tciclo = 1 / fmod;
nfit = 10; % Fit 10 cycles per volta
Tfit = nfit / fmod;

% Variables for all the fits
tmean = []; Acos = []; Asin = []; ADC = [];

% Medians in different time ranges
tf = []; Acosf = []; Asinf = [];
dAcosf = []; dAsinf = [];
ADCf = []; dADCf = [];

for kk = 1:size(tran, 1)
    N = floor((tran(kk, 2) - tran(kk, 1)) / Tfit);
    Acoskk = []; Asinkk = []; tmeankk = []; ADCkk = [];
    
    for nn = 1:N
        in = t > (nn-1) * Tfit + tran(kk, 1) & t <= nn * Tfit + tran(kk, 1);
        [fitout, dfitout] = fit_sine_poly(t(in), V(in), 0, fmod, 'center', 't0', tOFF(kk), 'nopl', 'nobs');
        
        ADCkk = [ADCkk; fitout(1)];
        Acoskk = [Acoskk; fitout(2)];
        Asinkk = [Asinkk; fitout(3)];
        tmeankk = [tmeankk; mean(t(in))];
    end
    
    Acos = [Acos; Acoskk];
    Asin = [Asin; Asinkk];
    ADC = [ADC; ADCkk];
    tmean = [tmean; tmeankk];
    
    Acosf = [Acosf; mean(Acoskk)];
    Asinf = [Asinf; mean(Asinkk)];
    dAcosf = [dAcosf; std(Acoskk) / sqrt(N-1)];
    dAsinf = [dAsinf; std(Asinkk) / sqrt(N-1)];
    ADCf = [ADCf; mean(ADCkk)];
    dADCf = [dADCf; std(ADCkk) / sqrt(N-1)];
    
    tf = [tf; mean(tmeankk)];
end

figdem = figure;
plot(tmean, Acos, 'b.');
hold on; grid on;
plot(tmean, Asin, 'r.');
legend('cos', 'sin', 'autoupdate', 'off');
xlabel('t (s)');
ylabel('V {DEMOD} (V)');

yl = [-0.6e-6 0.6e-6];
ylim(yl);

for jj = 1:size(tran, 1)
    plot(tran(jj, 1) * [1 1], yl, 'g-');
    plot(tran(jj, 2) * [1 1], yl, 'm--');
end

er = errorbar(tf, Acosf, dAcosf, 'cx');
set(er, 'linewidth', 3);

er = errorbar(tf, Asinf, dAsinf, 'mx');
set(er, 'linewidth', 3);
%}








