% TR = stlread('~/ecse/data/cube.stl');
% % this cube is centred on the origin and has side length 1
% % scaling is done first, then shifting
% % So for a regular array of cube size(s) size Hx, Hy, Hz
% % scale = [Hx Hy Hz], shift(n,:) = [2*n*Hx 2*n*Hy Hz/2];
% 
% ndivs = 4;
% TR_div = TR;
% for n=1:ndivs
%     TR_div = halfFaces(TR_div);
% end
% 
% 
% %%
% figure
% view(3)
% patch('Faces', TR_div.ConnectivityList, 'Vertices',TR_div.Points, 'FaceColor', ones(3,1)*69/100, 'FaceAlpha', 1)
% % hold on
% % incenters = TR_div.incenter;
% % faceNormals = TR_div.faceNormal;
% % quiver3(incenters(:,1), incenters(:,2), incenters(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3))
% axis equal

%%
function TR_out = divideFaces(TR_in)

conn = TR_in.ConnectivityList;
points = TR_in.Points;
npts = size(points,1);

conn_out = [];
pts_new = zeros(0,3);

npts_new = 0;

for n=1:size(conn,1)
    ptid1 = conn(n,1);
    ptid2 = conn(n,2);
    ptid3 = conn(n,3);
    pt1 = points(ptid1,:);
    pt2 = points(ptid2,:);
    pt3 = points(ptid3,:);

    % calculate distance between points
    dist12 = vecnorm(pt1 - pt2);
    dist13 = vecnorm(pt1 - pt3);
    dist23 = vecnorm(pt2 - pt3);
    dists = [dist23, dist13, dist12];
    [~,id] = max(dists);
    switch id
        case 1
            pt = (pt2 + pt3) / 2;
        case 2
            pt = (pt1 + pt3) / 2;
        case 3
            pt = (pt1 + pt2) / 2;
            
    end

    [lia, locb] = ismember(pt, pts_new, 'rows');
    if ~lia || isempty(pts_new)
        pts_new = [pts_new; pt];
        npts_new = npts_new + 1;
        ptid = npts + npts_new;

    else
        ptid = npts + locb;
    end

    switch id
        case 1
            conn_out = [conn_out; [ptid1 ptid2 ptid; ptid1 ptid ptid3]];
        case 2
            conn_out = [conn_out; [ptid1 ptid2 ptid; ptid2 ptid3 ptid]]; 
        case 3
            conn_out = [conn_out; [ptid1 ptid ptid3; ptid2 ptid3 ptid]]; 
             
    end 
end

pts_out = [points; pts_new];

TR_out = triangulation(conn_out, pts_out);

end


