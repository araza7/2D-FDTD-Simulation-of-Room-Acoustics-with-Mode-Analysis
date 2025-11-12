%% Script for 2D Finite Difference Time Domain (FDTD) simulation of room acoustics
% Author: Ali Raza

clear;      % Clear workspace variables
close all;  % Close all figures

%% --- Simulation Parameters ---

dt = 0.01 / 341 / sqrt(2); 
% Alternative: dt=0.01/341*(0.70); for a Courant number of 0.7

% Spatial step (dh) and grid dimensions
dh = 0.01;           % Spatial step in meters
nx = round(4 / dh);  % Number of grid points in x-direction (for a 4m room length)
ny = round(3 / dh);  % Number of grid points in y-direction (for a 3m room width)

% Medium Properties (Air)
rho = 1.21;         % Density of air (kg/m^3)
c = 341;            % Speed of sound in air (m/s)
k = (c^2) * rho;    % Bulk modulus (Pa)

% Absorption Coefficient and Impedance
alpha = 0.1;
R = sqrt(1 - alpha);         % Reflection coefficient magnitude
Znorm = (1 + R) / (1 - R);   % Normalized impedance for real R
impedance = Znorm * rho * c; % Specific acoustic impedance (Rayls)

%% --- Initialization of Acoustic Fields ---

p = zeros(nx, ny);       % Pressure field
ux = zeros(nx + 1, ny);  % Particle velocity in x-direction (staggered grid)
uy = zeros(nx, ny + 1);  % Particle velocity in y-direction (staggered grid)

%% --- Excitation Signal (Ricker Wavelet) ---

maxtt = 10000;               % Total number of time steps
centralfrequency = 1000;     % Central frequency of the Ricker wavelet (Hz)
a = centralfrequency / (sqrt(pi) / 2) * 4; % Ricker parameter

% Time vector and Ricker wavelet
t = ((1:maxtt) / (1/dt) - 4/a);
w = -(exp(-a^2 * (t.^2) / 2) .* (a^2 * (t.^2) - 1)); % Ricker wavelet

%% --- Data Recording ---

h2 = zeros(1, maxtt); % Preallocate array for output pressure at one point

%% --- Main Time-Stepping Loop ---

for tt = 1:maxtt
    % Update pressure field
    p = p - k * dt / dh * (diff(ux, 1, 1) + diff(uy, 1, 2));

    % Apply excitation (hard source) at p(1,1)
    p(1, 1) = w(tt);

    % Update velocity fields
    ux(2:nx, :) = ux(2:nx, :) - dt / rho / dh * diff(p, 1, 1);
    uy(:, 2:ny) = uy(:, 2:ny) - dt / rho / dh * diff(p, 1, 2);

    % Impedance boundary conditions (absorbing)
    ux(1,:)   = -p(1,:)/(rho*c*impedance);     % Left wall
    ux(end,:) =  p(end,:)/(rho*c*impedance);   % Right wall
    uy(:,1)   = -p(:,1)/(rho*c*impedance);     % Bottom wall
    uy(:,end) =  p(:,end)/(rho*c*impedance);   % Top wall

    % Record pressure at the opposite corner
    h2(tt) = p(nx, ny);
end

%% --- Post-processing: Frequency Analysis ---

% FFT of source and response
H1 = abs(fft(w, 100000));   % Input spectrum
H2 = abs(fft(h2, 100000));  % Output spectrum

% Frequency vector
n_fft = length(H1);
Fs = 1 / dt;
freq = (0:n_fft/2 - 1) * Fs / n_fft;

% Transfer function (output/input)
T_FDTD2 = H2 ./ (H1 + eps); % Avoid divide-by-zero
modes = 20 * log10(abs(T_FDTD2(1:n_fft/2)) + eps); % dB magnitude

%% --- Theoretical Room Modes Calculation ---

fprintf('Calculating theoretical room modes...\n');
lx = 4; % Room length (m)
ly = 3; % Room width (m)
ii = 0:100; jj = 0:100;
F = zeros(length(ii), length(jj));

for ii_loop = 0:100
    mx_modes = jj;
    my_mode = ii_loop;
    fc_calculated = c / 2 * sqrt(((mx_modes / lx).^2) + ((my_mode / ly).^2));
    F(ii_loop + 1, 1:length(jj)) = fc_calculated;
end

f = sort(F(:)'); % Flatten and sort all mode frequencies

%% --- Plotting Results ---

fprintf('Plotting results...\n');
figure;
semilogx(freq, modes, 'b', 'LineWidth', 2); hold on;

% Overlay theoretical mode lines
num_theoretical_modes_to_plot = 380;
plotted_modes_count = 0;

for kk = 1:length(f)
    if plotted_modes_count >= num_theoretical_modes_to_plot
        break;
    end
    if f(kk) > 1e-3
        if kk == 1 || abs(f(kk) - f(kk-1)) > 1e-3
            line([f(kk) f(kk)], [-3 20], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.4);
            plotted_modes_count = plotted_modes_count + 1;
        end
    end
end

% Plot formatting
xlim([20 1000]);
ylim([-24 20]);
xticks_values = [16 20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 ...
                 1250 1600 2000 2500 3150 4000 5000];
xticks(xticks_values);
xtickangle(45);

title('Comparison of Room Transfer Function and Theoretical Modes');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('FDTD Transfer Function', 'Theoretical Modes (Rigid Walls)');
grid on; hold off;
saveas(gcf, 'FDTD_Room_Transfer_Function.jpg');
