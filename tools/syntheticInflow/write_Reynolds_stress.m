%% Generate Input Profiles for Synthetic Inflow Turbulence Generator
% Two modes available:
%   Mode 1: Extract profiles from existing u-DALES simulation (NetCDF)
%   Mode 2: Generate custom profiles using empirical formulas
%
% Outputs: Reynolds stress profiles and length/time scale files

clear all; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ========== MODE SELECTION ==========

mode = 1;  % 1 = Extract from existing simulation, 2 = Custom profiles

%% ========== COMMON SETTINGS ==========

output_dir = 'E:\uDALES\jh6016c\u-dales-v2.02\u-dales\tools\syntheticInflow\driver_profiles';
stride = 4;          % Output every stride-th point
ltempeq = false;     % Include temperature profiles
lmoist = false;      % Include moisture profiles
lplot_profiles = true;

%% ========== MODE 1: EXTRACT FROM EXISTING SIMULATION ==========

if mode == 1
    
    % Source simulation settings
    expstr_source = '515';
    DA_SOURCE = 'D:/Postdoc1/simulation/testu2/outputs';
    iplane = 492;        % x-plane to extract from
    column1 = 4;         % y-range for averaging
    column2 = 7;
    
    % Read NetCDF file
    fileID = [DA_SOURCE '/' expstr_source '/tdump.' expstr_source '.nc'];
    
    zt = ncread(fileID, 'zt');
    zm = ncread(fileID, 'zm');
    ut = ncread(fileID, 'ut');
    upuptc = ncread(fileID, 'upuptc');
    vpvptc = ncread(fileID, 'vpvptc');
    wpwptc = ncread(fileID, 'wpwptc');
    upvpt = ncread(fileID, 'upvpt');
    upwpt = ncread(fileID, 'upwpt');
    vpwpt = ncread(fileID, 'vpwpt');
    
    nvec = length(zm);
    len = length(zm) + 1;
    
    z = zeros(len, 1);
    z(1:len-1) = zm;
    z(len) = zm(end) + (zm(end) - zm(end-1));
    
    % Extract profiles
    u = obtain_profile(ut, zt, iplane, column1, column2, z, zm, zt);
    upup = obtain_profile(upuptc, zt, iplane, column1, column2, z, zm, zt);
    vpvp = obtain_profile(vpvptc, zt, iplane, column1, column2, z, zm, zt);
    wpwp = obtain_profile(wpwptc, zt, iplane, column1, column2, z, zm, zt);
    upvp = obtain_profile(upvpt, zt, iplane, column1, column2, z, zm, zt);
    upwp = obtain_profile(upwpt, zm, iplane, column1, column2, z, zm, zt);
    vpwp = obtain_profile(vpwpt, zm, iplane, column1, column2, z, zm, zt);
    
    if ltempeq
        thlt = ncread(fileID, 'thlt');
        thlpthlpt = ncread(fileID, 'thlpthlpt');
        wpthlpt = ncread(fileID, 'wpthlpt');
        
        thl = obtain_profile(thlt, zt, iplane, column1, column2, z, zm, zt);
        thl(1) = thl(2);
        thlpthlp = obtain_profile(thlpthlpt, zt, iplane, column1, column2, z, zm, zt);
        wpthlp = obtain_profile(wpthlpt, zm, iplane, column1, column2, z, zm, zt);
        upthlp = wpthlp;
        vpthlp = wpthlp;
    end
    
    if lmoist
        qtt = ncread(fileID, 'qtt');
        qt = obtain_profile(qtt, zt, iplane, column1, column2, z, zm, zt);
        upqtp = wpthlp;
        vpqtp = upqtp;
        wpqtp = upqtp;
        thlpqtp = upqtp;
        qtpqtp = upqtp;
    end
    
    % Length and time scales (constant values)
    lengthy_u = 4 * ones(len, 1);
    lengthz_u = 4 * ones(len, 1);
    ts_u = 0.5 * ones(len, 1);
    lengthy_v = lengthy_u;
    lengthz_v = lengthz_u;
    ts_v = ts_u;
    lengthy_w = lengthy_u;
    lengthz_w = lengthz_u;
    ts_w = 0.5 * sqrt(2) * ones(len, 1);

%% ========== MODE 2: CUSTOM PROFILES (an example)==========

elseif mode == 2
    
    % Vertical grid
    z = (0:2:720)';
    len = length(z);
    
    % Mean velocity (power law profile)
    u_H = 11.7;      % Reference velocity [m/s]
    z_H = 180;       % Reference height [m]
    u = u_H * (z / z_H).^0.16;
    u(1) = 0;
    
    % Turbulence intensities
    I_H = 0.105;
    Iu = I_H * (z / z_H).^(-0.2);
    Iu(~isfinite(Iu)) = 0;
    Iv = 0.78 * Iu;
    Iw = 0.55 * Iu;
    
    % Reynolds stresses
    upup = (Iu .* u).^2;
    vpvp = (Iv .* u).^2;
    wpwp = (Iw .* u).^2;
    upvp = -0.3 * sqrt(upup .* vpvp);
    upwp = -0.2 * sqrt(upup .* wpwp);
    vpwp = -0.1 * sqrt(vpvp .* wpwp);
    
    % Length scales (grid points)
    dy = 4.0;
    dz_mean = mean(diff(z));
    L_H = 2.5 * z_H;
    Lu = L_H * (z / z_H).^0.133;
    Lv = 0.5 * Lu;
    Lw = 0.5 * Lu;
    
    lengthy_u = ceil(Lu / dy);
    lengthz_u = ceil(Lu / dz_mean);
    lengthy_v = ceil(Lv / dy);
    lengthz_v = ceil(Lv / dz_mean);
    lengthy_w = ceil(Lw / dy);
    lengthz_w = ceil(Lw / dz_mean);
    
    % Time scales
    ts_u = Lu ./ u;
    ts_v = Lv ./ u;
    ts_w = sqrt(ts_u.^2 + ts_v.^2);
    ts_u(~isfinite(ts_u)) = 0;
    ts_v(~isfinite(ts_v)) = 0;
    ts_w(~isfinite(ts_w)) = 0;
    
    % Fix surface values
    lengthy_u(1) = lengthy_u(2); lengthz_u(1) = lengthz_u(2); ts_u(1) = ts_u(2);
    lengthy_v(1) = lengthy_v(2); lengthz_v(1) = lengthz_v(2); ts_v(1) = ts_v(2);
    lengthy_w(1) = lengthy_w(2); lengthz_w(1) = lengthz_w(2); ts_w(1) = ts_w(2);
    
    % Optional: Temperature profiles
    if ltempeq
        zi = 1000;  % Boundary layer height
        thl = 300 + 5 * z / zi;
        thlpthlp = 0.1 * (1 - z/zi).^2;
        thlpthlp(z > zi) = 0.01;
        upthlp = 0.05 * sqrt(upup .* thlpthlp);
        vpthlp = 0.05 * sqrt(vpvp .* thlpthlp);
        wpthlp = -0.1 * sqrt(wpwp .* thlpthlp);
    end
    
    % Optional: Moisture profiles
    if lmoist
        zi = 1000;
        qt = 0.010 - 0.005 * z / zi;
        qt = max(qt, 0.001);
        qtpqtp = 1e-6 * (1 - z/zi).^2;
        upqtp = 0.01 * sqrt(upup .* qtpqtp);
        vpqtp = 0.01 * sqrt(vpvp .* qtpqtp);
        wpqtp = -0.05 * sqrt(wpwp .* qtpqtp);
        thlpqtp = ltempeq * (-0.5 * sqrt(thlpthlp .* qtpqtp));
    end
    
end

%% ========== CREATE OUTPUT DIRECTORY ==========

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
cd(output_dir);

%% ========== WRITE OUTPUT FILES ==========

% Velocity Reynolds stress
profile_vel = [z(1:stride:end) u(1:stride:end) upup(1:stride:end) upvp(1:stride:end) ...
               vpvp(1:stride:end) upwp(1:stride:end) vpwp(1:stride:end) wpwp(1:stride:end)];

filename = fopen('Reynolds_stress_profiles_velocity.txt', 'w');
fprintf(filename, "z u u'u' u'v' v'v' u'w' v'w' w'w'\n");
fprintf(filename, '%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', profile_vel');
fclose(filename);

fprintf('✓ Reynolds_stress_profiles_velocity.txt\n');

% Temperature (optional)
if ltempeq
    profile_temp = [z(1:stride:end) thl(1:stride:end) upthlp(1:stride:end) ...
                    vpthlp(1:stride:end) wpthlp(1:stride:end) thlpthlp(1:stride:end)];
    
    filename = fopen('Reynolds_stress_profiles_temp.txt', 'w');
    fprintf(filename, "z thl u'thl' v'thl' w'thl' thl'thl'\n");
    fprintf(filename, '%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', profile_temp');
    fclose(filename);
    
    fprintf('✓ Reynolds_stress_profiles_temp.txt\n');
end

% Moisture (optional)
if lmoist
    profile_moist = [z(1:stride:end) qt(1:stride:end) upqtp(1:stride:end) ...
                     vpqtp(1:stride:end) wpqtp(1:stride:end) thlpqtp(1:stride:end) qtpqtp(1:stride:end)];
    
    filename = fopen('Reynolds_stress_profiles_moist.txt', 'w');
    fprintf(filename, "z qt u'q' v'q' w'q' thl'q' q'q'\n");
    fprintf(filename, '%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n', profile_moist');
    fclose(filename);
    
    fprintf('✓ Reynolds_stress_profiles_moist.txt\n');
end

% Length and time scales
write_length_time_file('length_time_scales_u.txt', z, lengthy_u, lengthz_u, ts_u, stride);
write_length_time_file('length_time_scales_v.txt', z, lengthy_v, lengthz_v, ts_v, stride);
write_length_time_file('length_time_scales_w.txt', z, lengthy_w, lengthz_w, ts_w, stride);

fprintf('✓ length_time_scales_u.txt\n');
fprintf('✓ length_time_scales_v.txt\n');
fprintf('✓ length_time_scales_w.txt\n');

%% ========== PLOT PROFILES (OPTIONAL) ==========

if lplot_profiles
    profile_info = {'$z$ [m]', '$\overline{u}$', '$\overline{u^{\prime}u^{\prime}}$', ...
                    '$\overline{u^{\prime}v^{\prime}}$', '$\overline{v^{\prime}v^{\prime}}$', ...
                    '$\overline{u^{\prime}w^{\prime}}$', '$\overline{v^{\prime}w^{\prime}}$', ...
                    '$\overline{w^{\prime}w^{\prime}}$'};
    plot_profiles(profile_info, profile_vel);
    
    if ltempeq
        profile_info = {'$z$ [m]', '$\overline{\theta}$', '$\overline{u^{\prime}\theta^{\prime}}$', ...
                        '$\overline{v^{\prime}\theta^{\prime}}$', '$\overline{w^{\prime}\theta^{\prime}}$', ...
                        '$\overline{\theta^{\prime}\theta^{\prime}}$'};
        plot_profiles(profile_info, profile_temp);
    end
    
    if lmoist
        profile_info = {'$z$ [m]', '$\overline{q}$', '$\overline{u^{\prime}q^{\prime}}$', ...
                        '$\overline{v^{\prime}q^{\prime}}$', '$\overline{w^{\prime}q^{\prime}}$', ...
                        '$\overline{\theta^{\prime}q^{\prime}}$', '$\overline{q^{\prime}q^{\prime}}$'};
        plot_profiles(profile_info, profile_moist);
    end
end

fprintf('\nAll files written to: %s\n', output_dir);

%% ========== HELPER FUNCTIONS ==========

function profile = obtain_profile(var, ztzm, iplane, column1, column2, z, zm, zt)
    % Extract and interpolate profile from 4D NetCDF variable
    
    var_y_m = squeeze(mean(squeeze(mean(var(:, :, :, column1:column2), 4)), 2));
    vec = var_y_m(iplane, :)';
    
    nvec = length(zm);
    len = length(z);
    
    profile = zeros(len, 1);
    profile(2:len-1) = interp1(ztzm, vec, z(2:len-1), 'linear');
    profile(len) = vec(nvec) + ((vec(nvec) - vec(nvec-1)) / (ztzm(nvec) - ztzm(nvec-1))) * (z(len) - ztzm(nvec));
end

function write_length_time_file(filename, z, lengthy, lengthz, ts, stride)
    % Write length and time scale file
    
    fileID = fopen(filename, 'w');
    fprintf(fileID, 'z nLy nLz T\n');
    fprintf(fileID, '%15.10f %8d %8d %15.10f\n', ...
            [z(1:stride:end) lengthy(1:stride:end) lengthz(1:stride:end) ts(1:stride:end)]');
    fclose(fileID);
end

function plot_profiles(profile_info, profile)
    % Plot vertical profiles
    
    for i = 2:length(profile(1, :))
        figure;
        set(gca, 'Box', 'on', 'FontSize', 20.0, 'LineWidth', 2.5);
        grid on; hold on;
        plot(profile(:, i), profile(:, 1), 'LineWidth', 2);
        xlabel(profile_info{i}, 'FontSize', 30);
        ylabel(profile_info{1}, 'FontSize', 30);
        axis tight;
    end
end