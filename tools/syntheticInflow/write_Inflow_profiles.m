%% Custom Script to Write Synthetic Inflow Input Files
% This script allows manual input of Reynolds stress profiles and generates
% all necessary txt files for the synthetic inflow turbulence generator
%
% Required outputs:
%   - Reynolds_stress_profiles_velocity.txt
%   - Reynolds_stress_profiles_temp.txt (optional)
%   - Reynolds_stress_profiles_moisture.txt (optional)
%   - length_time_scales_u.txt
%   - length_time_scales_v.txt
%   - length_time_scales_w.txt

clear all; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ========== USER INPUT SECTION ==========

% Output directory
output_dir = 'E:\uDALES\jh6016c\work\pressure_flc\CAARC\770\syntheticInflow_inputs';

% Flags for optional outputs
ltempeq = false;  % Set to true if you want to generate temperature file
lmoist = false;   % Set to true if you want to generate moisture file

% Stride for output (output every stride-th point)
stride = 1;  % Set to 1 to output all points, 2 for every other point, etc.

%% ========== DEFINE VERTICAL LEVELS (z) ==========
% Option 1: Uniform spacing
% z = (0:2:200)';  % From 0 to 200m with 2m spacing

% Option 2: Non-uniform spacing (stretched grid)
% z_lower = (0:1:50)';     % 0-50m with 1m spacing
% z_upper = (55:5:200)';   % 50-200m with 5m spacing
% z = [z_lower; z_upper];

% Option 3: Manual input
z = [0:2:256]';
nz = length(z);

fprintf('Number of vertical levels: %d\n', nz);
fprintf('z range: %.2f to %.2f m\n', z(1), z(end));

%% ========== VELOCITY PROFILES ==========

% Mean velocity [m/s]
% Example: logarithmic profile
u_H = 6.0;       % Friction velocity [m/s]
z_H = 20;
%z0 = 0.1;           % Roughness length [m]
%kappa = 0.41;       % von Karman constant

% u = u_star / kappa * log((z + z0) / z0);
u = u_H * (z/ z_H).^0.16;
u(1) = 0;  % Set surface velocity to zero

figure
plot(u, z)
% OR manual input:
% u = [0; 1.5; 2.5; 3.2; 3.8; 4.3; 5.0; 5.5; 5.9; 6.2; 6.7; 7.0; 7.5; 7.8; 8.0; 8.1];

fprintf('\nVelocity statistics:\n');
fprintf('  u_mean = %.2f m/s\n', mean(u));
fprintf('  u_max  = %.2f m/s\n', max(u));

%% ========== REYNOLDS STRESS PROFILES ==========

% Reynolds stresses [m^2/s^2]
% Example: Empirical formulas based on boundary layer theory

% Friction velocity squared
% u_star2 = u_star^2;

% Initialize arrays
upup = zeros(nz, 1);
vpvp = zeros(nz, 1);
wpwp = zeros(nz, 1);
upvp = zeros(nz, 1);
upwp = zeros(nz, 1);
vpwp = zeros(nz, 1);

I_H = 0.105;
Iu = I_H*(z/z_H).^(-0.2);
Iu(~isfinite(Iu)) = 0;

Iv = 0.78* Iu;
Iw = 0.55* Iu;

upup = (Iu.*u).^2;
vpvp = (Iv.*u).^2;
wpwp = (Iw.*u).^2;
upvp = -0.3* sqrt(upup.*vpvp);
upwp = -0.2* sqrt(upup.*wpwp);
vpwp = -0.1* sqrt(vpvp.*wpwp);

fprintf('\nReynolds stress statistics:\n');
fprintf('  <u''u''>_max = %.4f m^2/s^2\n', max(upup));
fprintf('  <v''v''>_max = %.4f m^2/s^2\n', max(vpvp));
fprintf('  <w''w''>_max = %.4f m^2/s^2\n', max(wpwp));
fprintf('  <u''v''>_min = %.4f m^2/s^2\n', min(upvp));
figure
subplot(1,3,1); plot(Iu, z); xlabel('$I_u$')
subplot(1,3,2); plot(Iv, z); xlabel('$I_v$')
subplot(1,3,3); plot(Iw, z); xlabel('$I_w$')

% figure
% subplot(1,3,1); plot(upup, z/z_H); xlabel('$upup$')
% subplot(1,3,2); plot(vpvp, z/z_H); xlabel('$vpvp$')
% subplot(1,3,3); plot(wpwp, z/z_H); xlabel('$wpwp$')
% 
% figure
% subplot(1,3,1); plot(upvp, z/z_H); xlabel('$upvp$')
% subplot(1,3,2); plot(vpwp, z/z_H); xlabel('$vpwp$')
% subplot(1,3,3); plot(upwp, z/z_H); xlabel('$upwp$')

%% ========== LENGTH AND TIME SCALES ==========

% Grid spacing in target simulation (IMPORTANT: adjust these!)
dy = 2.0;  % y-direction grid spacing [m]
dz_mean = mean(diff(z));  % Average z-direction grid spacing [m]

L_H = 2.0 * z_H;
Lx_target = L_H * (z/z_H).^(0.133);
Ly_target =  0.3*Lx(floor(nz/2)) * ones(size(Lx_target));
Lz_target =  0.5*Lx(floor(nz/2)) * ones(size(Lx_target));

% Convert to real length, according to the filter constant
Lx = Lx_target*pi/2;
Ly = Ly_target*pi/2;
Lz = Lz_target*pi/2;

% Convert to grid points
lengthy_u = ceil(Ly/dy);           % Luuy = Ly
lengthz_u = ceil(Lz / dz_mean);    % Luuz = Lz
lengthy_v = ceil(Ly/dy);           % Lvvy = Ly
lengthz_v = ceil(Lz / dz_mean);    % Lvvz = Lz
lengthy_w = ceil(Ly/dy);           % Lwwy = Ly
lengthz_w = ceil(Lz / dz_mean);    % Lwwz = Lz

% Luux, Lvvx, Lwwx = Lx
ts_u = Lx ./ u;
ts_v = Lx ./ u;
ts_w = Lx ./ u;
ts_u(~isfinite(ts_u)) = 0;
ts_v(~isfinite(ts_v)) = 0;
ts_w(~isfinite(ts_w)) = 0;

lengthy_u(1) = lengthy_u(2); lengthz_u(1) = lengthz_u(2); ts_u(1) = ts_u(2);
lengthy_v(1) = lengthy_v(2); lengthz_v(1) = lengthz_v(2); ts_v(1) = ts_v(2);
lengthy_w(1) = lengthy_w(2); lengthz_w(1) = lengthz_w(2); ts_w(1) = ts_w(2);

figure
subplot(1,3,1); plot(Lx_target, z); xlabel('$L_x$')
subplot(1,3,2); plot(Ly_target, z); xlabel('$L_y$')
subplot(1,3,3); plot(Lz_target, z); xlabel('$L_z$')


%% ========== TEMPERATURE PROFILES (OPTIONAL) ==========

if ltempeq
    % Mean potential temperature [K]
    thl = 300 + 5 * z / zi;  % Linear profile with 5K across BL
    
    % Temperature variance [K^2]
    thlpthlp = 0.1 * (1 - z/zi).^2;
    thlpthlp(z > zi) = 0.01;
    
    % Temperature-velocity covariances [K*m/s]
    upthlp = 0.05 * sqrt(upup .* thlpthlp);
    vpthlp = 0.05 * sqrt(vpvp .* thlpthlp);
    wpthlp = -0.1 * sqrt(wpwp .* thlpthlp);  % Negative: upward heat flux
    
    fprintf('\nTemperature statistics:\n');
    fprintf('  thl: %.2f - %.2f K\n', min(thl), max(thl));
    fprintf('  <thl''thl''>_max = %.4f K^2\n', max(thlpthlp));
end

%% ========== MOISTURE PROFILES (OPTIONAL) ==========

if lmoist
    % Mean specific humidity [kg/kg]
    qt = 0.010 - 0.005 * z / zi;  % Decreasing with height
    qt = max(qt, 0.001);  % Minimum value
    
    % Moisture variance [(kg/kg)^2]
    qtpqtp = 1e-6 * (1 - z/zi).^2;
    
    % Moisture-velocity covariances [kg/kg * m/s]
    upqtp = 0.01 * sqrt(upup .* qtpqtp);
    vpqtp = 0.01 * sqrt(vpvp .* qtpqtp);
    wpqtp = -0.05 * sqrt(wpwp .* qtpqtp);
    
    % Temperature-moisture covariance [K * kg/kg]
    if ltempeq
        thlpqtp = -0.5 * sqrt(thlpthlp .* qtpqtp);
    else
        thlpqtp = zeros(nz, 1);
    end
    
    fprintf('\nMoisture statistics:\n');
    fprintf('  qt: %.5f - %.5f kg/kg\n', min(qt), max(qt));
    fprintf('  <qt''qt''>_max = %.2e (kg/kg)^2\n', max(qtpqtp));
end

%% ========== CREATE OUTPUT DIRECTORY ==========

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('\nCreated output directory: %s\n', output_dir);
end

cd(output_dir);

%% ========== WRITE VELOCITY REYNOLDS STRESS FILE ==========

profile_vel = [z(1:stride:end) u(1:stride:end) upup(1:stride:end) ...
               upvp(1:stride:end) vpvp(1:stride:end) upwp(1:stride:end) ...
               vpwp(1:stride:end) wpwp(1:stride:end)];

filename = fopen('Reynolds_stress_profiles_velocity.txt','w');
fprintf(filename,"z u u'u' u'v' v'v' u'w' v'w' w'w'\n");
fprintf(filename,'%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', profile_vel');
fclose(filename);

fprintf('\n✓ Written: Reynolds_stress_profiles_velocity.txt\n');

%% ========== WRITE TEMPERATURE FILE (IF ENABLED) ==========

if ltempeq
    profile_temp = [z(1:stride:end) thl(1:stride:end) upthlp(1:stride:end) ...
                    vpthlp(1:stride:end) wpthlp(1:stride:end) thlpthlp(1:stride:end)];
    
    filename = fopen('Reynolds_stress_profiles_temp.txt','w');
    fprintf(filename,"z thl u'thl' v'thl' w'thl' thl'thl'\n");
    fprintf(filename,'%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', profile_temp');
    fclose(filename);
    
    fprintf('✓ Written: Reynolds_stress_profiles_temp.txt\n');
end

%% ========== WRITE MOISTURE FILE (IF ENABLED) ==========

if lmoist
    profile_moist = [z(1:stride:end) qt(1:stride:end) upqtp(1:stride:end) ...
                     vpqtp(1:stride:end) wpqtp(1:stride:end) thlpqtp(1:stride:end) ...
                     qtpqtp(1:stride:end)];
    
    filename = fopen('Reynolds_stress_profiles_moisture.txt','w');
    fprintf(filename,"z qt u'qt' v'qt' w'qt' thl'qt' qt'qt'\n");
    fprintf(filename,'%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', profile_moist');
    fclose(filename);
    
    fprintf('✓ Written: Reynolds_stress_profiles_moisture.txt\n');
end

%% ========== WRITE LENGTH AND TIME SCALE FILES ==========

% For u-component
filename = 'length_time_scales_u.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'z nLy nLz T\n');
fprintf(fileID,'%15.10f %8d %8d %15.10f\n', ...
    [z(1:stride:end) lengthy_u(1:stride:end) lengthz_u(1:stride:end) ts_u(1:stride:end)]');
fclose(fileID);
fprintf('✓ Written: length_time_scales_u.txt\n');

% For v-component (same as u)
filename = 'length_time_scales_v.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'z nLy nLz T\n');
fprintf(fileID,'%15.10f %8d %8d %15.10f\n', ...
    [z(1:stride:end) lengthy_v(1:stride:end) lengthz_v(1:stride:end) ts_v(1:stride:end)]');
fclose(fileID);
fprintf('✓ Written: length_time_scales_v.txt\n');

% For w-component (same as u)
filename = 'length_time_scales_w.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'z nLy nLz T\n');
fprintf(fileID,'%15.10f %8d %8d %15.10f\n', ...
    [z(1:stride:end) lengthy_w(1:stride:end) lengthz_w(1:stride:end) ts_w(1:stride:end)]');
fclose(fileID);
fprintf('✓ Written: length_time_scales_w.txt\n');
