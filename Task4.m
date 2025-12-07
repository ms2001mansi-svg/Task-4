%% ========================================================================
% Project 3 – TASK 4 : Population PD Simulation 
%% ========================================================================

clear; clc; close all;

%% ============================
% INDEX ARRAY
%% ============================

n.VEGF_b = 1;  n.VEGFR2_b = 2;  n.VEGFVEGFR2_b = 3;  
n.Ab_b = 4;   n.VEGFAb_b = 5;

n.VEGF_t = 6;  n.VEGFR2_t = 7;  n.VEGFVEGFR2_t = 8;
n.Ab_t = 9;   n.VEGFAb_t = 10;

n.VEGF_r = 11; n.VEGFR2_r = 12; n.VEGFVEGFR2_r = 13;
n.Ab_r = 14;  n.VEGFAb_r = 15;

%% ============================
% PARAMETERS
%% ============================

p.Vol_b = 5 * 0.6;
p.Vol_t = 1 * 0.611;
p.Vol_r = 40 * 0.08;

p.kcl_V  = 0.0648;
p.kcl_VA = 2e-5;
p.kcl_A  = 2.2e-5;

p.kintR2  = 0.0168;     % name used in your ODE file
p.kintVR2 = p.kintR2;   % alias

p.VEGFR2prod_b = 0;
p.VEGFR2prod_t = 2.85e2 * p.kintR2;
p.VEGFR2prod_r = 2.2e3 * p.kintR2;

p.VEGFprod_b = 0;
p.VEGFprod_t = 10;
p.VEGFprod_r = 10;

p.konVR  = 6.4e4;
p.koffVR = 0.06;

p.konVA  = 5.52e-6;
p.koffVA = 0.012;

p.k_bt = 4.12e-4;
p.k_tb = 9.31e-4;
p.k_br = 1.66e-4;
p.k_rb = 3.44e-4;

p.AbEx = 1;

%% ============================
% INITIAL CONDITIONS
%% ============================

y0 = zeros(15,1);

%% ============================
% ODE15s — FAST SETTINGS
%% ============================

options = odeset('AbsTol',1e-6,'RelTol',1e-6,'MaxStep',2000);

%% ============================
% STEADY STATE
%% ============================
Tss = linspace(0, 60*24*10, 500);   % 10 days, 500 points ONLY → FAST
[T1, Y1] = ode15s(@VEGFAbeqns, Tss, y0, options, p, n);
ySS = Y1(end,:)';

%% ============================
% ANTIBODY DOSE
%% ============================

dose = 700e-3 / 150e3 * 1e12;   % 700 mg → mol → pM

yDose = ySS;
yDose(n.Ab_b) = yDose(n.Ab_b) + dose / p.Vol_b;

%% ============================
% FIXED TIME VECTOR FOR ALL PATIENTS 
%% ============================

endtime = 60*24*7;        % 7 days
Tspan = linspace(0, endtime, 600);   % only 600 points → very fast

%% ============================
% RUN ONE SIMULATION TO GET ntime
%% ============================

[Ttest, Ytest] = ode15s(@VEGFAbeqns, Tspan, yDose, options, p, n);
ntime = length(Ttest);

%% ============================
% LOAD VIRTUAL PATIENT FILE
%% ============================

data = readtable('VirtualPatientsTask4.csv');
VEGFprod_list   = data.sVEGF_tumor;
VEGFR2prod_list = data.sVEGFR2_Tumor;

numP = height(data);

%% ============================
% PREALLOCATE STORAGE
%% ============================

VEGF_t_all    = zeros(numP, ntime);
VEGFR2_t_all  = zeros(numP, ntime);
VEGFAb_t_all  = zeros(numP, ntime);

%% ============================
% POPULATION SIMULATION (FAST)
%% ============================

for i = 1:numP

    p.VEGFprod_t   = VEGFprod_list(i);
    p.VEGFR2prod_t = VEGFR2prod_list(i);

    yDose = ySS;
    yDose(n.Ab_b) = yDose(n.Ab_b) + dose / p.Vol_b;

    [Tp, Yp] = ode15s(@VEGFAbeqns, Tspan, yDose, options, p, n);

    VEGF_t_all(i,:)    = Yp(:, n.VEGF_t);
    VEGFR2_t_all(i,:)  = Yp(:, n.VEGFR2_t);
    VEGFAb_t_all(i,:)  = Yp(:, n.VEGFAb_t);
end

%% ============================
% PLOTS
%% ============================

figure;
plot(Tspan, VEGF_t_all', 'LineWidth', 1);
xlabel('Time (min)'); ylabel('VEGF (pM)');
title('Free VEGF in Tumor — 100 Virtual Patients');

figure;
plot(Tspan, VEGFR2_t_all', 'LineWidth', 1);
xlabel('Time (min)'); ylabel('VEGFR2 (pM)');
title('VEGFR2 in Tumor — 100 Virtual Patients');

figure;
plot(Tspan, VEGFAb_t_all', 'LineWidth', 1);
xlabel('Time (min)'); ylabel('VEGF-Ab Complex (pM)');
title('VEGF-Antibody Complex in Tumor — 100 Virtual Patients');

disp("SIMULATION COMPLETE — Fast mode using ode15s!");
