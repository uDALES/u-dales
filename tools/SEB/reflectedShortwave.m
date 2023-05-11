function Knet = reflectedShortwave(Sdir, Dsky, vf, svf, albedo)

Kin0 = (Sdir + Dsky * svf); % Available radiation

Knet0 = (1 - albedo) .* Kin0; % Initial net (absorbed) radiation
Kout0 = albedo .* Kin0; % Initial reflected radiation

Knet = Knet0; % Absorbed radiation - accumulated with iterations
Kout = Kout0; % Reflected radiation - updated each iteration

count = 0;
while true
    count = count + 1;
    vf_Kout = vf * Kout;
    Kadd = (1 - albedo) .* vf_Kout; % Additional absorbed radiation 
    Kout =      albedo  .* vf_Kout; % Updated reflected radiation
    Knet = Knet + Kadd;             % Updated net (absorbed) radiation 

    if (max(Kadd ./ (Knet - Kadd)) < 0.01)
        break
    end
end
