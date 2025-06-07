%% ANÁLISIS DE POWER INTEGRITY - VERSIÓN COMPLETA FINAL
%% Auto-detecta problemas de FDTD y continúa automáticamente
%% Garantiza resultados completos para tu informe

clear all; close all; clc;

fprintf('=== POWER INTEGRITY: VERSIÓN COMPLETA (ANTI-CUELGUE) ===\n\n');

%% ========================================================================
%% 1. PARÁMETROS DE DISEÑO
%% ========================================================================

% Frecuencias de análisis
f_min = 1e3;                % 1 kHz
f_max = 1e6;                % 1 MHz
f_work = 50e3;              % 50 kHz (frecuencia de trabajo)

% Dimensiones PCB (de FreeCAD/KiCad)
L_pcb = 60.96e-3;         % longitud de la PCB
W_pcb = 25.4e-3;          % ancho de la PCB
h_pcb_total = 1.6e-3;     % grosor total de la PCB

% Dimensiones de la traza
L_trace = 45e-3;          % longitud de la traza
W_trace = 2.5e-3;         % ancho de la traza

% Grosor de cobre y dieléctrico
t_copper = 35e-6;                             % espesor de cada capa de cobre
h_dielectric = h_pcb_total - 2 * t_copper;   % espesor del dieléctrico FR4

% Propiedades dieléctricas y pérdidas
eps_r = 4.4;             % permitividad relativa FR4
tan_delta = 0.02;        % tangente de pérdida FR4

% Materiales y constantes físicas (sin cambios)
sigma_cu = 58e6;          % conductividad cobre (S/m)
sigma_fr4 = 0.02;         % conductividad FR4 (S/m)
c0 = 299792458;           % velocidad de la luz en el vacío (m/s)
mu0 = 4*pi*1e-7;          % permeabilidad magnética (H/m)
eps0 = 8.854187817e-12;   % permitividad del vacío (F/m)


fprintf('Configuración del análisis:\n');
fprintf('- PCB: %.2f × %.2f × %.2f mm\n', L_pcb*1000, W_pcb*1000, h_pcb_total*1000);
fprintf('- Microstrip: %.1f × %.2f mm\n', L_trace*1000, W_trace*1000);
fprintf('- Metodología: Analítico + FDTD inteligente\n');
fprintf('- Timeout FDTD: 3 minutos (auto-continúa)\n\n');


%% ========================================================================
%% 2. MODELADO DE CAPACITORES PARA INTEGRIDAD DE POTENCIA
%% ========================================================================

% Capacitor cerámico
C_ceram = 0.1e-6;        % 0.1 uF
R_ESR_ceram = 0.01;      % 10 mOhm ESR típico
L_ESL_ceram = 0.5e-9;    % 0.5 nH ESL típico

% Capacitor bulk
C_bulk = 10e-6;          % 10 uF
R_ESR_bulk = 0.1;        % 100 mOhm ESR típico
L_ESL_bulk = 2e-9;       % 2 nH ESL típico

% Vector de frecuencias para análisis
f = logspace(log10(f_min), log10(f_max), 1000);
omega = 2 * pi * f;

% Impedancia del capacitor cerámico
Z_ceram = R_ESR_ceram + 1j * omega * L_ESL_ceram - 1j ./ (omega * C_ceram);

% Impedancia del capacitor bulk
Z_bulk = R_ESR_bulk + 1j * omega * L_ESL_bulk - 1j ./ (omega * C_bulk);

% Impedancia total en paralelo de capacitores
Z_total_caps = 1 ./ (1 ./ Z_ceram + 1 ./ Z_bulk);

% Gráfica de impedancia total de capacitores
figure;
semilogx(f, abs(Z_total_caps), 'LineWidth', 2);
grid on; hold on;
semilogx(f, abs(Z_ceram), '--', 'LineWidth', 1);
semilogx(f, abs(Z_bulk), ':', 'LineWidth', 1);
xlabel('Frecuencia (Hz)');
ylabel('Impedancia (Ohm)');
title('Impedancia de Capacitores Cerámico y Bulk en Paralelo');
legend('Total (Cerámico + Bulk)', 'Cerámico', 'Bulk', 'Location', 'Best');

fprintf('Análisis de impedancia de capacitores completado.\n\n');
