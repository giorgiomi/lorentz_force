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
clear; clc; close all;
addpath long_data\;

d = 1.115e-02; %in cm
VDC = [6.113 5.114 4.112 3.110 2.106 2.609 3.612];

Voltaggio_long=readmatrix("Voltaggio_lungo.csv");
Corrente_long=readmatrix("Corrente_lungo.csv");

Voltaggio_long=Voltaggio_long(:,2);
Corrente_long=Corrente_long(:,2);

%plot(1:length(Corrente_long),Corrente_long);

time=(0:(length(Voltaggio_long)-1))*0.5;

figure()
hold on;
grid on;
box on;

plot(time,Voltaggio_long,'LineWidth',0.2);

%%
figtimeseries = figure;
t=time';
V=Voltaggio_long;

plot(t, V); % Assuming it's a temporal series (t, V)
hold on; grid on;
vlim = [9.65e-3 9.7e-3];
ylim(vlim);

% Define intervals to analyze at different distances
tran = [840 11000; 12000 22000; 23000  34000; 35000 45000 ]; %15 h

figure(figtimeseries);
for jj = 1:size(tran, 1)
    plot(tran(jj, 1) * [1 1], vlim, 'g--');
    plot(tran(jj, 2) * [1 1], vlim, 'm--');
end

tOFF = [4; 0; 0; 0]; % Offset t0, manually set
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



