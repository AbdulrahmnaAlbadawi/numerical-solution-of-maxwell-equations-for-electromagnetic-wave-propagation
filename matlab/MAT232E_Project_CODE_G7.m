% --- 1D MAXWELL EQUATIONS: EXPLICIT FDTD VS. IMPLICIT CRANK-NICOLSON ---
% This script simulates a Gaussian pulse propagating in 1D free space.

clear; clc; clf;

% 1. PHYSICAL & COMPUTATIONAL PARAMETERS
c = 3e8;                    % Speed of light (m/s) [cite: 1, 128]
eps0 = 8.854e-12;           % Vacuum Permittivity [cite: 1, 19]
mu0 = 4*pi*1e-7;            % Vacuum Permeability [cite: 1, 16]
z_limit = 1.0;              % Length of the domain (meters)
dz = 0.005;                 % Spatial step (meters)
steps = round(z_limit/dz);  % Number of grid points [cite: 1, 71]
z = linspace(0, z_limit, steps);

% --- FDTD TIMING (Must satisfy CFL condition) ---
dt_fdtd = dz / c;           % Stable time step (CFL Limit) [cite: 1, 127]

% --- CRANK-NICOLSON TIMING (Breaking the CFL limit) ---
% Implicit methods are unconditionally stable, allowing larger steps [cite: 1, 200]
dt_cn = 4 * dt_fdtd;        % 4x larger time step than FDTD

% 2. INITIALIZE FIELDS
% FDTD Arrays (Staggered Yee Grid) [cite: 1, 77]
Ex_fdtd = zeros(1, steps);
Hy_fdtd = zeros(1, steps);

% Crank-Nicolson Arrays (Matrix based) [cite: 1, 162]
Ex_cn = zeros(steps, 1);
Ex_cn_old = zeros(steps, 1);

% 3. CRANK-NICOLSON MATRIX ASSEMBLY [cite: 1, 238]
% Solving the wave equation: E(n+1) - (c*dt)^2/2 * d2E/dz2(n+1) = RHS
s = (c * dt_cn / (2 * dz))^2;
main_diag = (1 + 2*s) * ones(steps, 1);
off_diag  = -s * ones(steps-1, 1);

% Matrix A is sparse and tridiagonal [cite: 1, 165]
A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

% Boundary Conditions: Perfect Electric Conductor (PEC) [cite: 1, 233]
A(1,:) = 0; A(1,1) = 1;
A(end,:) = 0; A(end,end) = 1;

% 4. MAIN SIMULATION LOOP
total_time = 150 * dt_fdtd;
n_steps = round(total_time / dt_fdtd);

for n = 1:n_steps
    current_time = n * dt_fdtd;

    % --- SIDE 1: EXPLICIT FDTD (Leapfrog) [cite: 1, 94] ---
    % Update H-field (Half-step) [cite: 1, 97]
    Hy_fdtd(1:end-1) = Hy_fdtd(1:end-1) + (dt_fdtd/(mu0*dz)) * diff(Ex_fdtd);
    
    % Update E-field (Full-step) [cite: 1, 99]
    % Corrected indexing: Matching lengths to avoid 'Incompatible Sizes' error
    Ex_fdtd(2:end-1) = Ex_fdtd(2:end-1) + (dt_fdtd/(eps0*dz)) * diff(Hy_fdtd(1:end-1));
    
    % Inject Gaussian Pulse Source [cite: 1, 230]
    pulse = exp(-((current_time - 20*dt_fdtd)/(6*dt_fdtd))^2);
    Ex_fdtd(round(steps/3)) = Ex_fdtd(round(steps/3)) + pulse;

    % --- SIDE 2: IMPLICIT CRANK-NICOLSON (Matrix Solve) [cite: 1, 235] ---
    % We update the CN solution only when its larger time step is reached
    if mod(n, 4) == 0
        % Calculate Right Hand Side (RHS) based on old values [cite: 1, 242]
        % For 1D wave: B * E_now
        rhs = 2*Ex_cn - Ex_cn_old + s*( [Ex_cn(2:end); 0] - 2*Ex_cn + [0; Ex_cn(1:end-1)] );
        
        % Inject Source into CN
        rhs(round(steps/3)) = rhs(round(steps/3)) + pulse;
        
        % Solve the system: A * E_next = RHS [cite: 1, 245]
        Ex_cn_next = A \ rhs; 
        
        % Shift time history
        Ex_cn_old = Ex_cn;
        Ex_cn = Ex_cn_next;
    end

    % --- VISUALIZATION ---
    if mod(n, 2) == 0
        plot(z, Ex_fdtd, 'b-', 'LineWidth', 1.5); hold on;
        plot(z, Ex_cn, 'r--', 'LineWidth', 1.5); hold off;
        grid on;
        axis([0 z_limit -1.5 2]);
        title(['Maxwell 1D: FDTD (Explicit) vs. Crank-Nicolson (Implicit)']);
        subtitle(['Time: ', num2str(current_time*1e9, '%.2f'), ' ns | Blue: FDTD | Red: CN (4x step)']);
        xlabel('Position (z)'); ylabel('Electric Field (Ex)');
        legend('FDTD (CFL Constrained)', 'Crank-Nicolson (Unconditionally Stable)');
        drawnow;
    end
end