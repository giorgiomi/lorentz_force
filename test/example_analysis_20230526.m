clear; close all; clc
filename = "Voltaggio_lungo.csv";

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 19);

% Specify range and delimiter
opts.DataLines = [10, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Index", "Reading", "Unit", "RangeDigits", "DispDigits", "Math", "StartGroup", "Limit1High", "Limit1Low", "Limit2High", "Limit2Low", "Terminal", "Questionable", "Origin", "Date", "Time", "FractionalSeconds", "Channel", "CHLabel"];
opts.VariableTypes = ["double", "double", "categorical", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "datetime", "datetime", "double", "categorical", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
% opts = setvaropts(opts, "CHLabel", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["Unit", "Math", "StartGroup", "Limit1High", "Limit1Low", "Limit2High", "Limit2Low", "Terminal", "Questionable", "Origin", "Channel", "CHLabel"], "EmptyFieldRule", "auto");
% opts = setvaropts(opts, "Date", "InputFormat", "MM/dd/yyyy");
% opts = setvaropts(opts, "Time", "InputFormat", "HH:mm:ss");

% ho dovuto commentare questo, forse non disponibile con il mio vecchio
% matlab (r2018b)

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
Index = tbl.Index;
V = tbl.Reading;
Unit = tbl.Unit;
Date = tbl.Date;
Time = tbl.Time;
FractionalSeconds = tbl.FractionalSeconds;

%% Clear temporary variables
clear opts tbl

%% Convert into seconds
utc = convertTo(Time, 'posixtime');
t = utc - utc(1) + FractionalSeconds;

%% somehow it did not manage the date correctly ... do this by hand ... 
day2_point1 = 33537;
day3_point1 = 212752;
n = [1:length(t)]'; 
day2 = n>=day2_point1; 
day3 = n>=day3_point1; 
t = t + day2*24*3600 + day3*24*3600; 


%% Prepare a plot
figure
plot(t, V)
hold on
grid on
xlabel("Time (s)")
ylabel(sprintf("V (%s)", Unit(1)))
%% 

fig_time_series = figure; 
plot(t,V);
hold on;
grid on;
vlim = [9.65e-3 9.7e-3];
ylim(vlim)
title("Data")

% definire intervalli da analizzare (ad esempio per diverse distanze)
t_ran = [1 11; 12 22; 23 34; 35 45;47 57;58 68;70 80;81 91;93 103;105 114;116 126;127 137; 139 200]*1e3;
fAC = [100; 200; 500; 1000;2000;3000;4000;50000;6000;7000;8000;10000;15000];
figure(fig_time_series);
for jj=1:size(t_ran,1)
    plot(t_ran(jj,1)*[1 1],vlim,'g--');
    plot(t_ran(jj,2)*[1 1],vlim,'m--');
end
%%
ind = [1:length(fAC)]';
t_OFF = [1.1*ones(size(t_ran,1))-ind*.00] ; % offset t_0 
% qui immagino un piccolo ritardo per l'impostazione oscillatore ogni ciclo ... 

f_mod = 0.1;            % frequenza modulazione forza
T_ciclo = 1/f_mod;          
n_fit = 10;             % analizzare fit a 10 cicli per volta
T_fit = n_fit/f_mod; 

% variabile per tutti i fit
t_mean = []; A_cos = []; A_sin = []; ADC = []; A_cos2 = []; A_sin2 = [];

% mediati sui gruppi nei diversi time range
t_f = []; Acos_f = []; Asin_f = []; dAcos_f = []; dAsin_f = [];
Acos2_f = []; Asin2_f = []; dAcos2_f = []; dAsin2_f = [];

ADC_f = []; dADC_f = [];

for kk = 1:size(t_ran,1)
    N = floor((t_ran(kk,2) - t_ran(kk,1))/ T_fit);
    A_coskk = [];    A_sinkk = [];    A_coskk2 = [];    A_sinkk2 = []; t_meankk = []; ADCkk = [];
    for nn = 1:N
        in = t>(nn-1)*T_fit + t_ran(kk,1) & t<=nn*T_fit + t_ran(kk,1);        
        [fit_out,dfit_out,C,chi2,N_DOF]=fit_sine_poly(t(in),V(in),0,f_mod*[1 2],'center','t0',t_OFF(kk),'nopl','nobs');
        ADCkk = [ADCkk; fit_out(1)]; 
        A_coskk = [A_coskk; fit_out(2)]; 
        A_sinkk = [A_sinkk; fit_out(3)];
        A_coskk2 = [A_coskk2; fit_out(4)]; 
        A_sinkk2 = [A_sinkk2; fit_out(5)];
        t_meankk = [t_meankk; mean(t(in))];
    end  
    A_cos = [A_cos; A_coskk];
    A_sin = [A_sin; A_sinkk];
    A_cos2 = [A_cos2; A_coskk2];
    A_sin2 = [A_sin2; A_sinkk2];
    ADC = [ADC; ADCkk];
    t_mean = [t_mean; t_meankk];
    
    Acos_f = [Acos_f; mean(A_coskk)];
    Asin_f = [Asin_f; mean(A_sinkk)];
    Acos2_f = [Acos2_f; mean(A_coskk2)];
    Asin2_f = [Asin2_f; mean(A_sinkk2)];
    dAcos_f = [dAcos_f; std(A_coskk)/sqrt(N-1)];
    dAsin_f = [dAsin_f; std(A_sinkk)/sqrt(N-1)];
    dAcos2_f = [dAcos2_f; std(A_coskk2)/sqrt(N-1)];
    dAsin2_f = [dAsin2_f; std(A_sinkk2)/sqrt(N-1)];
    ADC_f = [ADC_f; mean(ADCkk)];
    dADC_f = [dADC_f; std(ADCkk)/sqrt(N-1)];
    t_f = [t_f; mean(t_meankk)];
end

fig_dem = figure; 
plot(t_mean,A_cos,'b.');
hold on;
plot(t_mean,A_sin,'r.');
grid on; 
legend('cos','sin','autoupdate','off');
xlabel('t (s)');
ylabel('V_{DEMOD} (V)');
yl = [-0.6e-6 0.6e-6]; 
ylim(yl);
title('1f');

fig_dem2 = figure; 
plot(t_mean,A_cos2,'b.');
hold on;
plot(t_mean,A_sin2,'r.');
grid on; 
legend('cos','sin','autoupdate','off');
xlabel('t (s)');
ylabel('V_{DEMOD} (V)');
yl = [-0.6e-6 0.6e-6]; 
ylim(yl);
title('2f');

figure(fig_dem);
for jj=1:size(t_ran,1)
    plot(t_ran(jj,1)*[1 1], yl,'g-');
    plot(t_ran(jj,2)*[1 1], yl,'m--');
end
er = errorbar(t_f, Acos_f,dAcos_f,'cx'); set(er,'linewidth',3);
hold on
er = errorbar(t_f, Asin_f,dAsin_f,'mx'); set(er,'linewidth',3);


figure(fig_dem2);
for jj=1:size(t_ran,1)
    plot(t_ran(jj,1)*[1 1], yl,'g-');
    plot(t_ran(jj,2)*[1 1], yl,'m--');
end
er = errorbar(t_f, Acos2_f,dAcos2_f,'cx'); set(er,'linewidth',3);
er = errorbar(t_f, Asin2_f,dAsin2_f,'mx'); set(er,'linewidth',3);


%%
dV_dF = 2.4e-3; % usando valore trovato dai compagni per 12 V


figure ;
er = errorbar(fAC, Acos_f/dV_dF,dAcos_f/dV_dF,'gx'); set(er,'linewidth',3);
hold on
er = errorbar(fAC, Asin_f/dV_dF,dAsin_f/dV_dF,'mx'); set(er,'linewidth',3);
grid on; 
set(gca,'fontsize',16); 
set(gca,'xscale','log');
ylabel('F_{DEMOD} (N) ')
xlabel('Frequency (Hz)')
title('1f');

%return
% aggiungere modello (calcolato in un altro programma): 
plot(f_model,F_dem_model_1f_sin,'m');
legend('F_{1f} (cos)','F_{1f} (sin)','model sine 1f','autoupdate','off');
xlim([80 3500]);

figure ;
er = errorbar(fAC, Acos2_f/dV_dF,dAcos2_f/dV_dF,'gx'); set(er,'linewidth',3);
hold on
er = errorbar(fAC, Asin2_f/dV_dF,dAsin2_f/dV_dF,'mx'); set(er,'linewidth',3);
grid on; 
set(gca,'fontsize',16); 
set(gca,'xscale','log');
ylabel('F_{DEMOD} (N) ')
xlabel('Frequency (Hz)')
title('2f');
plot(f_model,F_dem_model_2f_cos,'g');
legend('F_{2f} (cos)','F_{2f} (sin)','model cosine 2f','autoupdate','off');
xlim([80 3500]);
