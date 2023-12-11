% This script is primarily designed to be used by write_inputs.m, but can
% also be called on its own to calculate shortwave.
% Assumes the following exist in the workspace
% fpath: where to write files to
% lwindows
% TR
% ltimedepsw
% xazimuth
% if ltimedepsw is true, then also need:
%   start (in datetime format), runtime, dtSP
%   longitude, latitude, timezone, elevation
% if ltimedepsw is false, then need lcustomSW
% if lcustomsw is false, then need start, longitude, latitude, timezone, elevation
% if lcustomsw is true, then need solarazimuth, solarzenith, irradiance, and Dsky
% if lscatter is true, the net shortwave is calculated using albedos, vf, svf, albedos


currentPath = pwd;
stk = dbstack; activeFilename = which(stk(1).file);
[folder, ~, ~] = fileparts(activeFilename);
addpath([folder '/SPA/'])

nfcts = size(TR.ConnectivityList,1);
npoints = size(TR.Points,1);

if ldirectShortwaveFortran
    % write geometry info
    fileID = fopen([fpath 'vertices.txt'],'w');
    fprintf(fileID,'%15.10f %15.10f %15.10f\n',TR.Points');
    fclose(fileID);
    fileID = fopen([fpath 'faces.txt'],'w');
    fprintf(fileID,'%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[TR.ConnectivityList TR.incenter TR.faceNormal]');
    fclose(fileID);

    cd(folder);
    system('gfortran -O3 directShortwave.f90 -o DS.exe');
    copyfile('DS.exe', fpath)
    delete DS.exe
    cd(currentPath)
end

if ~ltimedepsw
    if lcustomsw % custom solar position and irradiance - default true
        azimuth = solarazimuth - xazimuth;
    else
        sp = solarPosition(start, longitude, latitude, timezone, elevation);
        solarzenith = sp.zenith;
        azimuth = sp.azimuth - xazimuth;
        [ashraeA, ashraeB, ashraeC] = ASHRAE(start.Month);
        irradiance = ashraeA * exp(-ashraeB / cosd(solarzenith));
        Dsky = ashraeC * irradiance;
    end

    nsun = [sind(solarzenith)*cosd(azimuth), sind(solarzenith)*-sind(azimuth), cosd(solarzenith)];

    if ldirectShortwaveFortran
        disp('Calculating shortwave radiation using Fortran.')
        cd(fpath)
        writeInfo_directShortwave(nfcts, npoints, nsun, irradiance, resolution, fpath)
        if lwindows
            system('DS.exe');
        else
            system('./DS.exe');
        end
        Sdir = dlmread([fpath 'Sdir_fort.txt'], '', 0, 0);
        delete DS.exe;
        cd(currentPath)
    else
        disp('Calculating shortwave radiation using MATLAB.')
        show_plot_2d = false; % User-defined
        show_plot_3d = false;  % User-defined
        Sdir = directShortwave(F, V, nsun, irradiance, resolution, show_plot_2d, show_plot_3d);
        %dlmwrite('Sdir.txt', Sdir) % for debugging/visualisation purposes
        fileID = fopen([fpath 'Sdir.txt'], 'w');
        fprintf(fileID,'%8.2f\n', Sdir');
        fclose(fileID);
        disp('written Sdir.txt')
    end

    if lscatter
        % Calculate net shortwave radiation (Knet)
        Knet = netShortwave(Sdir, Dsky, vf, svf, albedos);
    end

else
    [ashraeA, ashraeB, ashraeC] = ASHRAE(start.Month);
    tSP = 0:dtSP:runtime;
    Sdir = zeros(nfcts, length(tSP));
    Knet = zeros(nfcts, length(tSP));
    for n = 1:length(tSP)
        n
        TOD = start + seconds(tSP(n));
        sp = solarPosition(TOD, longitude, latitude, timezone, elevation);
        solarzenith = sp.zenith;
        if solarzenith < 90
            azimuth = sp.azimuth - xazimuth;
            nsun = [sind(solarzenith)*cosd(azimuth), sind(solarzenith)*-sind(azimuth), cosd(solarzenith)];
            irradiance = ashraeA * exp(-ashraeB / cosd(solarzenith));
            Dsky = ashraeC * irradiance;
            if ldirectShortwaveFortran
                writeInfo_directShortwave(nfcts, npoints, nsun, irradiance, resolution, fpath)
                if lwindows
                    system('DS.exe');
                else
                    system('./DS.exe');
                end
                Sdir(:,n) = dlmread([fpath 'Sdir.txt'], '', 0, 0);
            else
                Sdir(:,n) = directShortwave(F, V, nsun, irradiance, 0.1, true, false);
            end

            if lscatter
                Knet(:,n) = reflectedShortwave(Sdir(:,n), Dsky, vf, svf, albedos);
            end
        end

    end
end


function writeInfo_directShortwave(nfcts, npoints, nsun, irradiance, resolution, fpath)
    fileID = fopen([fpath 'info_directShortwave.txt'],'w');
    fprintf(fileID,'%8d %8d\n',[nfcts, npoints]);
    fprintf(fileID,'%15.10f %15.10f %15.10f\n', nsun');
    fprintf(fileID,'%5.10f\n',irradiance);
    fprintf(fileID,'%5.10f\n',resolution);
    fclose(fileID);
end