function Knet = netShortwave(Sdir, Dsky, vf, svf, al)
% Calculates the net shortwave on a list of facets.
% Sdir: direct solar irradiance.
% Dsky: diffuse sky irradiance.
% vf: view factors.
% svf: sky view factors.
% al: albedo
% Knet: net shortwave irradiance.

Kin0 = (Sdir + Dsky * svf); % Available radiation

Knet0 = (1 - al) .* Kin0; % Initial net (absorbed) radiation
Kout0 = al       .* Kin0; % Initial reflected radiation

Knet = Knet0; % Absorbed radiation - accumulated with iterations
Kout = Kout0; % Reflected radiation - updated each iteration

count = 0;
while true
    count = count + 1;
    vf_Kout = vf * Kout;
    Kadd = (1 - al) .* vf_Kout; % Additional absorbed radiation
    Kout =      al  .* vf_Kout; % Updated reflected radiation
    Knet = Knet + Kadd;         % Updated net (absorbed) radiation

    if (max(Kadd ./ (Knet - Kadd)) < 0.01)
        break
    end
end
end
