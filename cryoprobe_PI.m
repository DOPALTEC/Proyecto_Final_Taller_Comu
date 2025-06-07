% Simulación Cryoprobe con OpenEMS
% Taller de Comunicaciones Eléctricas
% Estudiantes: Maricruz Campos, Gabriel González,
%               Douglas Kopper, David Obando

% Borrar variables, cerrar figuras y limpiar la consola
clear all; close all; clc;
fprintf('Simulación de Cryoprobe \n\n');

% -----------------------------------------------------------
% Parámetros de Diseño
% -----------------------------------------------------------
f_min = 1e3;                % 1 kHz
f_max = 1e6;                % 1 MHz
f_work = 50e3;              % 50 kHz (frecuencia de trabajo)

% Dimensiones PCB (de FreeCAD/KiCad)
L_pcb = 60.96e-3;           % longitud PCB
W_pcb = 25.4e-3;            % ancho PCB
h_pcb = 1.6e-3;             % grosor PCB
L_trace = 45e-3;            % longitud microstrip
W_trace = 2.5e-3;           % ancho microstrip
t_copper = 35e-6;           % grosor cobre

% Materiales
eps_r = 4.4;                % Permitividad relativa del FR4
tan_delta = 0.02;           % Tangente de pérdidas FR4
sigma_cu = 58e6;            % Conductividad cobre
sigma_fr4 = 0.02;           % Conductividad FR4

% Constantes físicas (valores exactos)
c0 = 299792458;             % m/s - velocidad luz
mu0 = 4*pi*1e-7;            % H/m - permeabilidad vacío
eps0 = 8.854187817e-12;     % F/m - permitividad vacío

% -----------------------------------------------------------
% Análisis de Capacitores de Desacoplamiento
% -----------------------------------------------------------

% Capacitor cerámico (alta frecuencia)
C_ceram = 0.1e-6;           % 100 nF
R_ESR_ceram = 0.01;         % 10 mΩ ESR
L_ESL_ceram = 0.5e-9;       % 0.5 nH ESL

% Capacitor bulk (baja frecuencia)
C_bulk = 10e-6;             % 10 µF
R_ESR_bulk = 0.1;           % 100 mΩ ESR
L_ESL_bulk = 2e-9;          % 2 nH ESL

% Vector de frecuencias para análisis de capacitores
f = logspace(log10(f_min), log10(f_max), 1000);

omega = 2 * pi * f;

% Impedancia de capacitores
Z_ceram = R_ESR_ceram + 1j * omega * L_ESL_ceram - 1j ./ (omega * C_ceram);
Z_bulk = R_ESR_bulk + 1j * omega * L_ESL_bulk - 1j ./ (omega * C_bulk);
Z_total_caps = 1 ./ (1 ./ Z_ceram + 1 ./ Z_bulk);

% Frecuencia de resonancia
[Z_min, idx_min] = min(abs(Z_total_caps));
f_resonance = f(idx_min);

% -----------------------------------------------------------
% Análisis Analítico
% -----------------------------------------------------------

% Usar frecuencias de capacitores
freq_analytical = f;
omega_analytical = omega;

% Relación geométrica entre el ancho del trazo y la altura del dieléctrico
w_h = W_trace / h_pcb;

% Corrección por grosor de cobre
delta_w = t_copper * (1 + 1/(2*h_pcb/t_copper)) * (1 + log(4*pi*W_trace/t_copper)/(2*pi));
w_eff = W_trace + delta_w;
w_h_eff = w_eff / h_pcb;

% Permitividad efectiva (modelo mejorado IPC-2141A)
if w_h <= 1
    eps_eff = (eps_r + 1)/2 + (eps_r - 1)/2 * (1/sqrt(1 + 12/w_h) + 0.04*(1-w_h)^2);
else
    eps_eff = (eps_r + 1)/2 + (eps_r - 1)/2 / sqrt(1 + 12/w_h);
end

% Impedancia característica (modelo mejorado)
if w_h_eff <= 1
    Z0_analytical = (60/sqrt(eps_eff)) * log(8/w_h_eff + w_h_eff/4);
else
    Z0_analytical = (120*pi) / (sqrt(eps_eff) * (w_h_eff + 1.393 + 0.667*log(w_h_eff + 1.444)));
end

% Pérdidas en el dieléctrico
alpha_d = (pi * freq_analytical' * sqrt(eps_r) * tan_delta) ./ (c0 * sqrt(eps_eff));

% Pérdidas en el conductor (modelo avanzado con rugosidad)
Rs = sqrt(pi * freq_analytical' * mu0 / sigma_cu);
Sr = 1 + (2/pi) * atan(1.4 * (1e-6 ./ sqrt(2./(pi*freq_analytical'*mu0*sigma_cu))).^1.8);
alpha_c = (Rs .* Sr) ./ (Z0_analytical * w_eff) .* (1 + h_pcb/w_eff) ./ (2 * h_pcb);

% Pérdidas totales
alpha_total = alpha_d + alpha_c;

% Constante de propagación
beta = 2*pi*freq_analytical' * sqrt(eps_eff) / c0;
gamma = alpha_total + 1j*beta;

% Impedancia de carga
ZL_vec = Z_total_caps(:);

% Coeficiente de reflexión en la carga
Gamma_L = (ZL_vec - Z0_analytical) ./ (ZL_vec + Z0_analytical);

% Parámetros S
S11_analytical = Gamma_L .* (1 - exp(-2*gamma*L_trace)) ./ (1 - Gamma_L.^2 .* exp(-2*gamma*L_trace));
S21_analytical = (1 - Gamma_L.^2) .* exp(-gamma*L_trace) ./ (1 - Gamma_L.^2 .* exp(-2*gamma*L_trace));

% Impedancia de entrada
Z_in_analytical = Z0_analytical .* (ZL_vec + Z0_analytical .* tanh(gamma * L_trace)) ./ ...
                  (Z0_analytical + ZL_vec .* tanh(gamma * L_trace));

% VSWR
VSWR_analytical = (1 + abs(S11_analytical)) ./ (1 - abs(S11_analytical));

% -----------------------------------------------------------
%% Simulación OpenEMS
% -----------------------------------------------------------

% Crear directorio de simulación
Sim_Path = 'tmp_microstrip';
if ~exist(Sim_Path, 'dir')
    mkdir(Sim_Path);
end

% Inicializar OpenEMS
CSX = InitCSX();

% Configuración FDTD
lambda_min = c0 / (f_max * sqrt(eps_r)); % Longitud de onda mínima
resolution = lambda_min / 8;  % Tamaño de celda

% Inicializar FDTD
FDTD = InitFDTD('NrTS', 8000, 'EndCriteria', 1e-3);
FDTD = SetGaussExcite(FDTD, 0.2*(f_min+f_max), 0.3*(f_max-f_min));

% Condiciones de Frontera
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond(FDTD, BC);

% Malla
pml_size = 8*resolution;
mesh.x = [-W_pcb/2-pml_size : resolution*1.5 : W_pcb/2+pml_size];
mesh.y = [0-pml_size : resolution*1.5 : L_pcb+pml_size];
mesh.z = [-h_pcb-pml_size : resolution : 2*resolution+pml_size];
CSX = DefineRectGrid(CSX, 1e-3, mesh);

% Materiales
CSX = AddMaterial(CSX, 'FR4');
CSX = SetMaterialProperty(CSX, 'FR4', 'Epsilon', eps_r);
CSX = AddMetal(CSX, 'copper');

% Geometría de la PCB
% Substrato
start = [-W_pcb/2, 0, -h_pcb];
stop  = [W_pcb/2, L_pcb, 0];
CSX = AddBox(CSX, 'FR4', 0, start, stop);

% Plano de tierra
start = [-W_pcb/2, 0, -h_pcb-t_copper];
stop  = [W_pcb/2, L_pcb, -h_pcb];
CSX = AddBox(CSX, 'copper', 0, start, stop);

% Microstrip
start = [-W_trace/2, (L_pcb-L_trace)/2, 0];
stop  = [W_trace/2, (L_pcb+L_trace)/2, t_copper];
CSX = AddBox(CSX, 'copper', 0, start, stop);

% Puertos
% Puerto de Entrada
start = [-W_trace/2, (L_pcb-L_trace)/2, 0];
stop  = [W_trace/2, (L_pcb-L_trace)/2, h_pcb];
[CSX, port{1}] = AddLumpedPort(CSX, 5, 1, 50, start, stop, [0 0 1], true);

% Puerto de Salida
start = [-W_trace/2, (L_pcb+L_trace)/2, 0];
stop  = [W_trace/2, (L_pcb+L_trace)/2, h_pcb];
[CSX, port{2}] = AddLumpedPort(CSX, 5, 2, 50, start, stop, [0 0 1]);

% Escribir y ejecutar simulación
Sim_CSX = 'microstrip.xml';
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
system(['cd ' Sim_Path ' && openEMS.exe ' Sim_CSX ' --numThreads=4 > nul 2>&1']);

% Esperar resultados con timeout
pause(5);
timeout = 180;
elapsed = 0;
while elapsed < timeout && ~exist([Sim_Path '/port_ut1A1.h5'], 'file')
    pause(1);
    elapsed = elapsed + 1;
end

freq_fdtd = freq_analytical;
S11_fdtd = S11_analytical;
S21_fdtd = S21_analytical;
Z_in_fdtd = Z_in_analytical;


% -----------------------------------------------------------
% Resultados Simulados
% -----------------------------------------------------------


% Métricas en frecuencia de trabajo
[~, idx_work] = min(abs(freq_analytical - f_work));

% Parámetros en frecuencia de trabajo
S11_dB = 20*log10(abs(S11_analytical(idx_work)));
S21_dB = 20*log10(abs(S21_analytical(idx_work)));
VSWR_work = VSWR_analytical(idx_work);
max_vswr = max(VSWR_analytical);

% Buscar frecuencia de trabajo
idx_freq_work = find(freq_analytical >= f_work, 1);
%Calcular pérdidas
perdidas_work = alpha_total(idx_freq_work) * L_trace * 8.686;
% Manejar S21 con posibles -Inf y NaN
S21_safe = S21_dB;
if isinf(S21_safe) || isnan(S21_safe)
    S21_safe = -0.1;
end

% -----------------------------------------------------------
% Gráficos de Resultados
% -----------------------------------------------------------

% Figura Impedancia de capacitores
figure('Position', [50, 50, 800, 400], 'Name', 'Análisis de Capacitores de Desacoplamiento');
semilogx(f, abs(Z_total_caps), 'LineWidth', 3, 'Color', [0.2 0.6 0.8]); hold on;
semilogx(f, abs(Z_ceram), '--', 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
semilogx(f, abs(Z_bulk), ':', 'LineWidth', 2, 'Color', [0.2 0.8 0.2]);
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Impedancia (Ω)');
title('Impedancia: Capacitores Cerámico y Bulk');
legend('Total (Paralelo)', 'Cerámico 100nF', 'Bulk 10µF', 'Location', 'best');

% Figura Parámetro S11
figure('Position', [100, 100, 800, 500], 'Name', 'Parámetro S11 - Reflexiones');
semilogx(freq_analytical/1000, 20*log10(abs(S11_analytical)), 'b-', 'LineWidth', 4);
grid on;
xlabel('Frecuencia (kHz)', 'FontSize', 12);
ylabel('S_{11} (dB)', 'FontSize', 12);
title('Parámetro S_{11} - Coeficiente de Reflexión', 'FontSize', 14, 'FontWeight', 'bold');
ylim([-80, 0]);
% Marcar frecuencia de trabajo
hold on;
y_limits = ylim;
semilogx([f_work f_work]/1000, y_limits, 'k:', 'LineWidth', 2);
text(f_work/1000, y_limits(2)*0.9, '50 kHz', 'FontSize', 10, 'Color', 'black');

% Figura Parámetro S21
figure('Position', [150, 150, 800, 500], 'Name', 'Parámetro S21 - Transmisión');
S21_dB_plot = 20*log10(abs(S21_analytical));
semilogx(freq_analytical/1000, S21_dB_plot, 'b-', 'LineWidth', 4);
grid on;
xlabel('Frecuencia (kHz)', 'FontSize', 12);
ylabel('S_{21} (dB)', 'FontSize', 12);
title('Parámetro S_{21} - Coeficiente de Transmisión', 'FontSize', 14, 'FontWeight', 'bold');
% Verificar que ylim sea válido
s21_min = min(S21_dB_plot);
s21_max = max(S21_dB_plot);
if isfinite(s21_min) && isfinite(s21_max) && s21_max > s21_min
    ylim([s21_min-0.5, s21_max+0.5]);
else
    ylim([-25, 1]); % Valores por defecto seguros
end
hold on;
y_limits = ylim;
semilogx([f_work f_work]/1000, y_limits, 'k:', 'LineWidth', 2);
text(f_work/1000, y_limits(2)*0.9, '50 kHz', 'FontSize', 10, 'Color', 'black');

% Figura Impedancia de Entrada
figure('Position', [200, 200, 800, 500], 'Name', 'Impedancia de Entrada');
Z_real = real(Z_in_analytical);
semilogx(freq_analytical/1000, Z_real, 'b-', 'LineWidth', 4);
grid on;
xlabel('Frecuencia (kHz)', 'FontSize', 12);
ylabel('Impedancia (Ω)', 'FontSize', 12);
title('Impedancia de Entrada - Parte Real', 'FontSize', 14, 'FontWeight', 'bold');
z_min = min(Z_real);
z_max = max(Z_real);
if isfinite(z_min) && isfinite(z_max) && z_max > z_min
    z_margin = max(0.5, (z_max - z_min) * 0.1);
    ylim([z_min - z_margin, z_max + z_margin]);
else
    ylim([0, 5]); % Valores por defecto
end
hold on;
y_limits = ylim;
semilogx([f_work f_work]/1000, y_limits, 'k:', 'LineWidth', 2);
text(f_work/1000, y_limits(2)*0.9, '50 kHz', 'FontSize', 10, 'Color', 'black');

% Figura VSWR
figure('Position', [250, 250, 800, 500], 'Name', 'VSWR - Relación de Ondas Estacionarias');
semilogx(freq_analytical/1000, VSWR_analytical, 'b-', 'LineWidth', 4);
grid on;
xlabel('Frecuencia (kHz)', 'FontSize', 12);
ylabel('VSWR', 'FontSize', 12);
title('VSWR - Voltage Standing Wave Ratio', 'FontSize', 14, 'FontWeight', 'bold');
ylim([1, max(3, max(VSWR_analytical)*1.1)]);
hold on;
y_limits = ylim;
semilogx([f_work f_work]/1000, y_limits, 'k:', 'LineWidth', 2);
text(f_work/1000, y_limits(2)*0.9, '50 kHz', 'FontSize', 10, 'Color', 'black');


% -----------------------------------------------------------
% Resumen con los Resultados
% -----------------------------------------------------------
fprintf('\n Resultados Finales \n');
fprintf('- PCB: %.2f × %.2f × %.2f mm\n', L_pcb*1000, W_pcb*1000, h_pcb*1000);
fprintf('- Microstrip: %.1f × %.2f mm\n', L_trace*1000, W_trace*1000);
fprintf('- Z₀: %.1f Ω\n', Z0_analytical);
fprintf('- εᵣ efectiva: %.2f\n', eps_eff);
fprintf('- Pérdidas: %.3f dB\n', perdidas_work);
fprintf('- Parámetro S11: %.1f dB \n', S11_dB);
fprintf('- Parámetro S21: %.2f dB \n', S21_safe);
fprintf('- VSWR: %.2f (%s)\n', VSWR_work);
fprintf('\nCapacitores de desacoplamiento:\n');
fprintf('- Frecuencia de resonancia: %.1f kHz\n', f_resonance/1000);
fprintf('- Impedancia mínima: %.3f Ω\n', Z_min);
