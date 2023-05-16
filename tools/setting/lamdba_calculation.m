%% Lambda Calculator

Area = X*Y;

con = TR.ConnectivityList;
ps = TR.Points;
norms = TR.faceNormal;
west_area = 0;
east_area = 0;
green_area = 0;
top_area = 0;
green_facets = [];
for i = 1:length(con(:,1))
    verts =con(i,:);
    p1 = ps(verts(1),:);
    p2 = ps(verts(2),:);
    p3 = ps(verts(3),:);
    zsum = p1(3) + p2(3) + p3(3);
    if zsum == 0
        for g = 1:length(green_outline(:,1))
            gsp = green_outline(g,:);
            gxl = gsp(1);
            gxu = gsp(2);
            gyl = gsp(3);
            gyu = gsp(4);
            if all([and(p1(1)>=gxl,p1(1)<=gxu),...
                and(p1(2)>=gyl,p1(2)<=gyu),...
                and(p2(1)>=gxl,p2(1)<=gxu),...
                and(p2(2)>=gyl,p2(2)<=gyu),...
                and(p3(1)>=gxl,p3(1)<=gxu),...
                and(p3(2)>=gyl,p3(2)<=gyu)])
                area = area_facets(i);
                green_area = green_area + area;
                green_factets = [green_facets; i];
            end
        end
    else 
        fnorm = norms(i,:);
        if fnorm == [-1,0,0]
            area = area_facets(i);
            west_area = west_area + area;

        elseif fnorm == [1,0,0]
            area = area_facets(i);
            east_area = east_area + area;
        elseif fnorm == [0,0,1]
            area = area_facets(i);
            top_area = top_area + area;
        end     

    end
end
lp = top_area/Area;
lv = green_area/Area;
lf = west_area/Area;
