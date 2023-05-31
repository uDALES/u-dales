function max_facet_side = max_facet_size(vertices,facets)

% fpath = 'D:/Postdoc1/simulation/testu2/experiments/501/';
% TR = stlread([fpath 'ground.stl']);
% vertices = TR.Points;
% facets = TR.ConnectivityList;
% incenters = TR.incenter;
% faceNormals = TR.faceNormal;

max_facet_side = 0;
side = [];

for i=1:length(facets(:,1))
    side(1) = sqrt( (vertices(facets(i,2),1)-vertices(facets(i,1),1))^2 ...
                  + (vertices(facets(i,2),2)-vertices(facets(i,1),2))^2 ...
                  + (vertices(facets(i,2),3)-vertices(facets(i,1),3))^2 );
    
    side(2) = sqrt( (vertices(facets(i,3),1)-vertices(facets(i,2),1))^2 ...
                  + (vertices(facets(i,3),2)-vertices(facets(i,2),2))^2 ...
                  + (vertices(facets(i,3),3)-vertices(facets(i,2),3))^2 );
              
    side(3) = sqrt( (vertices(facets(i,1),1)-vertices(facets(i,3),1))^2 ...
                  + (vertices(facets(i,1),2)-vertices(facets(i,3),2))^2 ...
                  + (vertices(facets(i,1),3)-vertices(facets(i,3),3))^2 );
    
    if (max_facet_side<max(side))  
        max_facet_side = max(side);
    end
end