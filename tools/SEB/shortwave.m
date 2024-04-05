% This script is primarily designed to be used by write_inputs.m, but can
% also be called on its own to calculate shortwave.
% Assumes the following exist in the workspace
% fpath: where to write files to
% lwindows
% TR
% ltimedepsw
% xazimuth
% if ltimedepsw is true, then need:
% start (in datetime format), runtime, dtSP
% if lweatherfile is true; weatherfname
% else need longitude, latitude, timezone, elevation
% if ltimedepsw is false, can also use lcustomsw
% if lcustomsw is true, then need solarazimuth, solarzenith, irradiance,
%   and Dsky (and not longitude etc or a weatherfile)
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
end

if ~lscatter
    % geometric factor accounting for facet orientation for diffuse sky irradiance.
    Fss = (1 + dot(TR.faceNormal, repmat([0 0 1], nfcts, 1), 2)) / 2;
end

cd(fpath)
if ~ltimedepsw
    if lcustomsw % custom solar position and irradiance
        azimuth = solarazimuth - xazimuth;
    elseif lweatherfile
        datestring = datestr(start, 'ddmmyy');
        timenum = start.Hour * 3600;
        timedepdata = readtable(weatherfname);
        id = find(table2array(timedepdata(:,1)) == str2num(datestring) & ...
            table2array(timedepdata(:,2)) == timenum);
        Tair = table2array(timedepdata(id, 'TAIR')) + 273.15; % not used here, but in case needed
        irradiance = table2array(timedepdata(id, 'HELIOM'));
        Dsky = table2array(timedepdata(id, 'DIFSOLAR'));
        Lsky = table2array(timedepdata(id, 'IRSKYT')); % not used here, but in case needed
        solarzenith = table2array(timedepdata(id, 'SOLAR'));
        solarazimuth = table2array(timedepdata(id, 'SOLAR_1'))+90;
        azimuth = solarazimuth - xazimuth;
    else
        sp = solarPosition(start, longitude, latitude, timezone, elevation);
        solarzenith = sp.zenith;
        azimuth = sp.azimuth - xazimuth;
        [ashraeA, ashraeB, ashraeC] = ASHRAE(start.Month);
        irradiance = ashraeA * exp(-ashraeB / cosd(solarzenith));
        Dsky = ashraeC * irradiance;
    end

    if (irradiance > 0 || Dsky > 0)
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
            Sdir = dlmread([fpath 'Sdir.txt'], '', 0, 0);
            delete vertices.txt faces.txt info_directShortwave.txt DS.exe;
            %cd(currentPath)
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
        else
            Knet = (1 - albedos) .* (Sdir + Dsky * Fss);
        end

    else
        Sdir = zeros(nfcts, 1);
        Knet = zeros(nfcts, 1);
        if ldirectShortwaveFortran
            delete vertices.txt faces.txt info_directShortwave.txt DS.exe;
        end
    end

else
    if lweatherfile
        datestring = datestr(start, 'ddmmyy');
        timedepdata = readtable(weatherfname);
        ids = find(table2array(timedepdata(:,1)) == str2num(datestring));
        timedeptime = table2array(timedepdata(ids, 'TIME'));
        timedepTair = table2array(timedepdata(ids, 'TAIR')) + 273.15;
        timedepI = table2array(timedepdata(ids, 'HELIOM'));
        timedepDsky = table2array(timedepdata(ids, 'DIFSOLAR'));
        timedepLsky = table2array(timedepdata(ids, 'IRSKYT'));
        timedepzenith = table2array(timedepdata(ids, 'SOLAR'));
        timedepazimuth = table2array(timedepdata(ids, 'SOLAR_1'))+90;
        
        timedepazimuth_shift = circshift(timedepazimuth, -hour(start));
        timedepzenith_shift = circshift(timedepzenith, -hour(start));
        timedepI_shift = circshift(timedepI, -hour(start));
        timedepDsky_shift = circshift(timedepDsky, -hour(start));

        tSP = 0:dtSP:runtime;
        % currently assumes the weather file is cyclic
        zenith_interp = interp1([timedeptime; 86400], [timedepzenith_shift; timedepzenith_shift(1)], tSP, 'makima');
        azimuth_interp = interp1([timedeptime; 86400], [timedepazimuth_shift; timedepazimuth_shift(1)], tSP, 'makima');
        I_interp = interp1([timedeptime; 86400], [timedepI_shift; timedepI_shift(1)], tSP, 'makima');
        Dsky_interp = interp1([timedeptime; 86400], [timedepDsky_shift; timedepDsky_shift(1)], tSP, 'makima');

        %%
        Sdir = zeros(nfcts, length(tSP));
        Knet = zeros(nfcts, length(tSP));
        for n = 1:length(tSP)
            n
            solarzenith = zenith_interp(n);
            irradiance = I_interp(n);
            if (solarzenith < 90 && irradiance > 0)
                solarazimuth = azimuth_interp(n);
                azimuth = solarazimuth - xazimuth;
                nsun = [sind(solarzenith)*cosd(azimuth), sind(solarzenith)*-sind(azimuth), cosd(solarzenith)];
                Dsky = Dsky_interp(n);
                if ldirectShortwaveFortran
                    writeInfo_directShortwave(nfcts, npoints, nsun, irradiance, resolution, fpath)
                    if lwindows
                        system('DS.exe');
                    else
                        system('./DS.exe');
                    end
                    Sdir(:,n) = dlmread([fpath 'Sdir.txt'], '', 0, 0);
                else
                    Sdir(:,n) = directShortwave(F, V, nsun, irradiance, resolution, true, false);
                end

                if lscatter
                    Knet(:,n) = netShortwave(Sdir(:,n), Dsky, vf, svf, albedos);
                end
            end
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
                    Sdir(:,n) = directShortwave(F, V, nsun, irradiance, resolution, false, false);
                end

                if lscatter
                    Knet(:,n) = netShortwave(Sdir(:,n), Dsky, vf, svf, albedos);
                else
                    Knet(:,n) = (1 - albedos) .* (Sdir(:,n) + Dsky * Fss);
                end
            end

        end
    end
    if ldirectShortwaveFortran
        delete vertices.txt faces.txt info_directShortwave.txt DS.exe Sdir.txt;
    end

    %% write to netcdf for visualisation
    ncid = netcdf.create(['Sdir.nc'], 'NC_WRITE');
    dimidrow = netcdf.defDim(ncid,'rows', nfcts);
    dimidcol = netcdf.defDim(ncid,'columns', length(tSP));
    varid_tSP = netcdf.defVar(ncid,'tSP','NC_FLOAT',dimidcol);
    varid_Sdir = netcdf.defVar(ncid,'Sdir','NC_FLOAT',[dimidrow dimidcol]);
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid_tSP,tSP);
    netcdf.putVar(ncid,varid_Sdir,Sdir);
    netcdf.close(ncid);
end
cd(currentPath)

function writeInfo_directShortwave(nfcts, npoints, nsun, irradiance, resolution, fpath)
    fileID = fopen([fpath 'info_directShortwave.txt'],'w');
    fprintf(fileID,'%8d %8d\n',[nfcts, npoints]);
    fprintf(fileID,'%15.10f %15.10f %15.10f\n', nsun');
    fprintf(fileID,'%5.10f\n',irradiance);
    fprintf(fileID,'%5.10f\n',resolution);
    fclose(fileID);
end
