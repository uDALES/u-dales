function A = facetAreas(F, V)

nF = size(F,1);
A = zeros(nF,1);

for i = 1:nF
    verts = V(F(i,:), :);
    nml = cross(verts(2,:)-verts(1,:),verts(3,:)-verts(1,:));
    A(i) = 0.5 * norm(nml);
end
end