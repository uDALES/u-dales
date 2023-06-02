function S = directShortwave(F, V, nsun, irradiance, resolution, show_plot_2d, show_plot_3d)
% Calculates the direct shortwave radiation S using polygon scan
%   conversion.
% F (#facets x 3/4) is the connectivity list of facets.
% V (#vertices x 3) is the corresponding vertices.
% Together these describe the geometry, which should not contain any facets 
%   that are not included in the energy balance (apart from bounding walls).
%   Note the geometry can be composed of triangles or quadrilaterals.
% nsun (1 x 3) is the vector TO the sun (doesn't have to be unitary).
% irradiance [W/m^2] is the radiation on an object whose surface normal is
%   in the direction of nsun.

% [nv, Nf] = size(F);
[Nf, nv] = size(F);

if nv == 3
    TR = triangulation(F,V);
    % Can use the included methods of the class.
elseif nv == 4
    
else
    error('Only triangles or quadrilaterals allowed.')
end

%% Find centres
if nv == 3
    C = TR.incenter;
elseif nv == 4
    C = zeros(Nf, 3);
    for i = 1:Nf
        C(i,:) = sum(V(F(i,:), :), 1) / nv;
    end
end



%% Find surface normals
if nv == 3
    N = TR.faceNormal;
elseif nv == 4
    N = zeros(Nf,3);
    for i = 1:Nf
        % Assumes anticlockwise vertex ordering
        t1 = V(F(1,i),:) - V(F(2,i),:);
        t2 = V(F(1,i),:) - V(F(3,i),:);
        n = cross(t1,t2);
        N(i,:) = n / norm(n);
    end
end



%% Find areas
A = zeros(Nf, 1);

if nv == 3
    for i = 1:size(TR.ConnectivityList,1)
        verts = V(F(i,:), :);
        nml = cross(verts(2,:)-verts(1,:),verts(3,:)-verts(1,:));
        A(i) = 0.5 * norm(nml);
    end
elseif nv == 4
    for i = 1:Nf
        % Assumes vertex ordering
        A(i) = norm(V(F(i,1),:) - V(F(i,2),:)) * norm(V(F(i,1),:) - V(F(i,4),:));
    end
end

%% Projection
nsun = nsun / norm(nsun);

% Self-shading logical vector
% If SS == true then it is NOT self shading
SS = N * nsun' > 0;

% Find centre of geometry at the ground
Lx = max(V(:,1)) - min(V(:,1));
Ly = max(V(:,2)) - min(V(:,2));
p0 = [min(V(:,1)) + Lx/2, min(V(:,2)) + Ly/2, 0] + 3 * max(Lx, Ly) * nsun;

% Find orthonormal basis on plane
d = dot(nsun, p0); % Perpendicular distance of plane from origin
u1 = [nsun(2) -nsun(1) 0]; % Horizontal basis vector
u1 = -u1;
u1 = u1 / norm(u1);
u2 = cross(u1, nsun); % Vertical basis vector

%% Project
M = [u1; u2; -nsun];

% Centres
% Store local coordinates xi, eta in Cpl
% Store global coordinates x, y, z in Cpg
% Store distances to plane in mu

sC = (C - p0) / M;
Cpl  = sC(:,1:2);
muC  = sC(:,3);
Cpg = C + muC * nsun;

sV = (V - p0) / M;
Vpl  = sV(:,1:2);
muV  = sV(:,3);
Vpg = V + muV * nsun; 

mu = zeros(Nf,1);
for i=1:Nf
    %v = F(i,:);
    %mu(i) = min([muV(v), muC(i)]);
    mu(i) = muC(i);
end

% Sort in descending order
[~, I] = sort(mu, 1, 'descend');
ISS = I(SS(I) == 1);

% Shift local coordinates so that they are all positive
Vplshift = Vpl + abs(min(Vpl,[],1));

delta = resolution;

res_xi = ceil(max(Vplshift(:,1)) / delta);
res_eta = ceil(max(Vplshift(:,2)) / delta);
bw = zeros(res_eta, res_xi);

%figure

n=0;
for i = I'
    n=n+1;
    disp(['Surface: ' num2str(n) ' ; ~ ' num2str(round(n/Nf * 100, 1)) ' % complete'])
    mask = poly2mask(Vplshift(F(i,:), 1) / delta, Vplshift(F(i,:), 2) / delta, res_eta, res_xi);
    if SS(i) == 0
        bw(mask) = 0;
    else
        bw(mask) = mask(mask) * i;
    end
%   bw(mask) = mask(mask) * i * SS(i);

end

%% Find radiation
Ap = histc(bw(:), 1:Nf) * delta^2;
S = irradiance * Ap ./ A; % Radiation on facet (W/m2)

if show_plot_2d  
    Sbm = zeros(size(bw));
    for i = I'
        Sbm(bw == i) = S(i);
        %Sbm(bw == i) = i;
    end

    figure
    imagesc(Sbm)
    colormap gray
    c = colorbar;
    c.Title.String = 'Direct solar radiation [W/m^2]';
    %caxis([0,1])
    xlabel('\xi')
    ylabel('\eta')
    axis equal tight
end
    
if show_plot_3d
    figure
    grid on
%     xlim([0 plotlim])
%     ylim([0 plotlim])
%     zlim([0 plotlim])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    
    % Plot in descending order of distance to plane in order to plot
    % projection onto plane too
    %for i = I
    for i = I'
        clr =  repmat(S(i), 1, 3) / max(S);
        %clr = [1,1,1];
        %patch(V(1, F(:,i)), V(2, F(:,i)), V(3, F(:,i)), clr)
        patch(V(F(i,:),1), V(F(i,:),2), V(F(i,:),3), clr)
        %scatter3(C(1,i), C(2,i), C(3,i))
        %plot3([C(1,i), C(1,i) + N(1,i)],  [C(2,i), C(2,i) + N(2,i)],  [C(3,i), C(3,i) + N(3,i)]);
        %patch(Vpg(1, F(:,i)), Vpg(2, F(:,i)), Vpg(3, F(:,i)), clr, 'LineStyle', 'none')
        patch(Vpg(F(i,:),1), Vpg(F(i,:),2), Vpg(F(i,:),3), clr, 'LineStyle', 'none')
    end
%      plot3([origin(1) origin(1) + nsun(1)],  [origin(2) origin(2) + nsun(2)],  [origin(3) origin(3) + nsun(3)]);
%      plot3([origin(1) origin(1) + u1(1)], [origin(2) origin(2) + u1(2)], [origin(3) origin(3) + u1(3)]);
%      plot3([origin(1) origin(1) + u2(1)],[origin(2) origin(2) + u2(2)],[origin(3) origin(3) + u2(3)]);
     hold off
     axis equal tight
     view(nsun)
end
%end


