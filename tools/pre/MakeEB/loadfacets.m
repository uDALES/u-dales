function [fct, wall] = loadfacets(expnr)
M = dlmread(['facets.inp.' num2str(expnr)],'',1,0);
vars = {'or', 'wlid', 'blk', 'bld'};
for n = 1:length(vars)
    fct.(vars{n}) = M(:, n);
end

M = dlmread(['walltypes.inp.' num2str(expnr)],'',3,0);
% wallid  lGR   z0 [m]     z0h [m]     al [-]   em [-]   d1 [m]  d2 [m]    cp1 [J/kgK]  cp2 [J/kgK]   l1 [W/(m K)]  l2 [W/(m K)]    k1 [W/(m K)]    k1 [W/(m K)]
vars = {'id','lGR','z0','z0h','al', 'em', 'd1', 'd2', 'cp1', 'cp2', 'l1', 'l2', 'k1', 'k2'};
for n = 1:length(vars)
    wall.(vars{n}) = M(:, n);
end

% check whether all referenced wall ids have been defined
idref = unique(fct.wlid);
iddef = unique(wall.id');

if length(wall.id) > length(iddef)
    disp('ERROR: multiple definitions of a wall id')
    return
end
if length(idref) > length(iddef)
    disp('ERROR: more walltypes used than defined')
    return
end


% assign the wall type to each of the facets -- this way it is easy to know
% in which structure the proporties can be found.
fct.wltp = zeros(size(fct.or));

for i=1:size(fct.wlid,1)
    if fct.wlid(i) >= 0 && fct.wlid(i)<=10 && any(ismember(iddef,fct.wlid(i)))
        fct.wltp(n)=1; %normal
    elseif fct.wlid(i) > 10 && any(ismember(iddef,fct.wlid(i)))
        fct.wltp(n)=2; %Green roof
    elseif fct.wlid(i)<0 && any(ismember(iddef,fct.wlid(i)))
        if fct.wlid(i)<=-99
            fct.wltp(n)=3; %bounding wall
        else
            fct.wltp(n)=4; %floor
        end
    else
        disp(['ERROR: walltype of facet ',num2str(i),' not defined'])
    end
end


for n = 1:length(idref)
    if (ismember(idref(n), wall.id))
        fct.wltp(fct.wlid == idref(n)) = 1;
    elseif (ismember(idref(n), gr.id))
        fct.wltp(fct.wlid == idref(n)) = 2;
    else
        sprintf('ERROR. wall id %8d not defined. Correct wall and greenroof property input files', idref(n))
        return
    end
end