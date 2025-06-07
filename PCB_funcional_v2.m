%% ANÁLISIS DE POWER INTEGRITY - VERSIÓN COMPLETA FINAL (LIMPIO)
clear all; close all; clc;
fprintf('=== POWER INTEGRITY: VERSIÓN COMPLETA ===\n\n');

%% 1. PARÁMETROS DE DISEÑO
f_min = 1e3;
f_max = 1e6;
f_work = 50e3;

L_pcb = 60.96e-3;
W_pcb = 25.4e-3;
h_pcb_total = 1.6e-3;

L_trace = 45e-3;
W_trace = 2.5e-3;

t_copper = 35e-6;
h_dielectric = h_pcb_total - 2 * t_copper;

eps_r = 4.4;
tan_delta = 0.02;

sigma_cu = 58e6;
sigma_fr4 = 0.02;
c0 = 299792458;
mu0 = 4*pi*1e-7;
eps0 = 8.854187817e-12;

%% 2. MODELO DE CAPACITORES
C_ceram = 0.1e-6;
R_ESR_ceram = 0.01;
L_ESL_ceram = 0.5e-9;

C_bulk = 10e-6;
R_ESR_bulk = 0.1;
L_ESL_bulk = 2e-9;

f = logspace(log10(f_min), log10(f_max), 1000);
omega = 2 * pi * f;

Z_ceram = R_ESR_ceram + 1j * omega * L_ESL_ceram - 1j ./ (omega * C_ceram);
Z_bulk = R_ESR_bulk + 1j * omega * L_ESL_bulk - 1j ./ (omega * C_bulk);
Z_total_caps = 1 ./ (1 ./ Z_ceram + 1 ./ Z_bulk);

figure;
semilogx(f, abs(Z_total_caps), 'LineWidth', 2); hold on;
semilogx(f, abs(Z_ceram), '--', 'LineWidth', 1);
semilogx(f, abs(Z_bulk), ':', 'LineWidth', 1);
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Impedancia (Ohm)');
title('Impedancia: Capacitores Cerámico y Bulk');
legend('Total', 'Cerámico', 'Bulk');

fprintf('Análisis de capacitores completado.\n\n');

%% PARTE I: ANÁLISIS ANALÍTICO (IPC-2141A)
fprintf('PARTE I: ANÁLISIS ANALÍTICO\n');
tic_analytical = tic;

freq_analytical = f;
omega_analytical = omega;
w_h = W_trace / h_pcb_total;

if w_h <= 1
    eps_eff = (eps_r + 1)/2 + (eps_r - 1)/2 * (1/sqrt(1 + 12/w_h) + 0.04*(1-w_h)^2);
else
    eps_eff = (eps_r + 1)/2 + (eps_r - 1)/2 / sqrt(1 + 12/w_h);
end

delta_w = t_copper * (1 + 1/(2*h_pcb_total/t_copper)) * (1 + log(4*pi*W_trace/t_copper)/(2*pi));
w_eff = W_trace + delta_w;
w_h_eff = w_eff / h_pcb_total;

if w_h_eff <= 1
    Z0_analytical = (60/sqrt(eps_eff)) * log(8/w_h_eff + w_h_eff/4);
else
    Z0_analytical = (120*pi) / (sqrt(eps_eff) * (w_h_eff + 1.393 + 0.667*log(w_h_eff + 1.444)));
end

alpha_d = (pi * freq_analytical' * sqrt(eps_r) * tan_delta) ./ (c0 * sqrt(eps_eff));
Rs = sqrt(pi * freq_analytical' * mu0 / sigma_cu);
Sr = 1 + (2/pi) * atan(1.4 * (1e-6 ./ sqrt(2./(pi*freq_analytical'*mu0*sigma_cu))).^1.8);
alpha_c = (Rs .* Sr) ./ (Z0_analytical * w_eff) .* (1 + h_pcb_total/w_eff) ./ (2 * h_pcb_total);
alpha_total = alpha_d + alpha_c;

beta = 2*pi*freq_analytical' * sqrt(eps_eff) / c0;
gamma = alpha_total + 1j*beta;

ZL_vec = Z_total_caps(:); % Impedancia de carga en función de frecuencia
Gamma_L = (ZL_vec - Z0_analytical) ./ (ZL_vec + Z0_analytical);

S11_analytical = Gamma_L .* (1 - exp(-2*gamma*L_trace)) ./ (1 - Gamma_L.^2 .* exp(-2*gamma*L_trace));
S21_analytical = (1 - Gamma_L.^2) .* exp(-gamma*L_trace) ./ (1 - Gamma_L.^2 .* exp(-2*gamma*L_trace));

Z_in_analytical = Z0_analytical .* (ZL_vec + Z0_analytical .* tanh(gamma * L_trace)) ./ ...
                  (Z0_analytical + ZL_vec .* tanh(gamma * L_trace));
VSWR_analytical = (1 + abs(S11_analytical)) ./ (1 - abs(S11_analytical));

time_analytical = toc(tic_analytical);
fprintf('Tiempo: %.3f s\n', time_analytical);
fprintf('Z0: %.1f Ω\n', Z0_analytical);
fprintf('ε_eff: %.2f\n', eps_eff);

[~, idx_work_analytical] = min(abs(freq_analytical - f_work));
fprintf('Z_in @ 50kHz: %.1f Ω\n', abs(Z_in_analytical(idx_work_analytical)));
fprintf('VSWR @ 50kHz: %.2f\n', VSWR_analytical(idx_work_analytical));
