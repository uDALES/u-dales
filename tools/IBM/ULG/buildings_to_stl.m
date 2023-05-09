xs = [0:X];
ys = [0:Y];

points = [];

for i = 1:length(xs)
    for j = 1:length(ys)
        points =[points; xs(i), ys(j)];
    end 
end
for i = length(points(:,1)):-1:1
    p = points(i,:);
    px = p(1);
    py = p(2);
    for b = 1:length(buildings(:,1))
        bxl = buildings(b,1);
        bxu = buildings(b,2);
        byl = buildings(b,3);
        byu = buildings(b,4);
        bzu = buildings(b,6);
        if bzu > 0
            if and(px>bxl, px<bxu)
                if and(py>byl, py<byu)
                    points(i,:) = [];
                end 
            end
        end      
    end 
end

conlist = delaunay(points)


for b = 1:length(buildings(:,1))
    bxl = buildings(b,1);
    bxu = buildings(b,2);
    byl = buildings(b,3);
    byu = buildings(b,4);
    bzu = buildings(b,6);
    bplanx = [bxl,bxu];
    bplany = [byl,byu];
    if bzu> 0
        for c = length(conlist(:,1)):-1:1
            lp1 = 0;
            lp2 = 0;
            lp3 = 0;
            verts = conlist(c,:);
            p1 = points(verts(1),:);
            p2 = points(verts(2),:);
            p3 = points(verts(3),:);
            if and(p1(1)<=bxu,p1(1)>=bxl)
                if and(p1(2)<=byu,p1(2)>=byl)
                    disp('p1 lies on edge of building')
                    lp1 = 1;
                end 
            end
            if and(p2(1)<=bxu,p2(1)>=bxl)
                if and(p2(2)<=byu,p2(2)>=byl)
                    disp('p2 lies on edge of building')
                    lp2 = 1;
                end 
            end
            if and(p3(1)<=bxu,p3(1)>=bxl)
                if and(p3(2)<=byu,p3(2)>=byl)
                    disp('p3 lies on edge of building')
                    lp3 = 1;
                end 
            end
            prod = lp1.*lp2.*lp3;
            if prod == 1
                disp('all points on edge so delete.')
                conlist(c,:) = [];
            end 
        end 
    end     
end         
figure()
triplot(conlist, points(:,1),points(:,2))
xsf = points(:,1);
ysf = points(:,2);
zsf =zeros(size(points(:,1)))
conf = conlist;

%% Make building triangulation

xs = [];
ys = [];
zs = [];
con = [];
tot = length(xsf);
for b = 1:length(buildings(:,1))
    block = buildings(b,:)
    xl = block(1);
    xu = block(2);
    yl = block(3);
    yu = block(4);
    zl = block(5);
    zu = block(6);
    if zu ~= 0;
        xsb = xl:1:xu;
        ysb = yl:1:yu;
        zsb = zl:1:zu;
        %% Top
        pointstop = [];
        for i = 1:length(xsb)
            for j = 1:length(ysb)
                pointstop =[pointstop; xsb(i), ysb(j)];
            end 
        end
        contop = delaunay(pointstop);
        contop = contop + tot;
        tot = tot + length(pointstop(:,1));
        conf = [conf;contop];
        xsf = [xsf;pointstop(:,1)];
        ysf = [ysf;pointstop(:,2)];
        zsf = [zsf;zu.*ones(size(pointstop(:,1)))]; 

        %% West
        pointsw = [];
        for i = 1:length(zsb)
            for j = 1:length(ysb)
                pointsw =[pointsw; zsb(i), ysb(j)];
            end 
        end
        conw = delaunay(pointsw);
        conw = conw + tot;
        tot = tot + length(pointsw(:,1));
        conf = [conf;conw];
        xsf = [xsf;xl.*ones(size(pointsw(:,1)))];
        ysf = [ysf;pointsw(:,2)];
        zsf = [zsf;pointsw(:,1)];

        %% East
        pointse = [];
        for i = 1:length(zsb)
            for j = 1:length(ysb)
                pointse =[pointse; zsb(i), ysb(j)];
            end 
        end
        cone = delaunay(pointse);
        contemp = cone;
        cone(:,2) = contemp(:,3);
        cone(:,3) = contemp(:,2);
        cone = cone + tot;
        tot = tot + length(pointse(:,1));
        conf = [conf;cone];
        xsf = [xsf;xu.*ones(size(pointse(:,1)))];
        ysf = [ysf;pointse(:,2)];
        zsf = [zsf;pointse(:,1)];

        %% North
        pointsn = [];
        for i = 1:length(xsb)
            for j = 1:length(zsb)
                pointsn =[pointsn; xsb(i), zsb(j)];
            end 
        end
        conn = delaunay(pointsn);
        contemp = conn;
        conn(:,2) = contemp(:,3);
        conn(:,3) = contemp(:,2);
        conn = conn + tot;
        tot = tot + length(pointsn(:,1));
        conf = [conf;conn];
        xsf = [xsf;pointsn(:,1)];
        ysf = [ysf;yu.*ones(size(pointsn(:,1)))];
        zsf = [zsf;pointsn(:,2)];

        %% South
        pointss = [];
        for i = 1:length(xsb)
            for j = 1:length(zsb)
                pointss =[pointss; xsb(i), zsb(j)];
            end 
        end
        cons = delaunay(pointss);
        cons = cons + tot;
        tot = tot + length(pointss(:,1));
        conf = [conf;cons];
        xsf = [xsf;pointss(:,1)];
        ysf = [ysf;yl.*ones(size(pointss(:,1)))];
        zsf = [zsf;pointss(:,2)];

    end
end 
TR = triangulation(conf,xsf,ysf,zsf);
%% Plot
figure
%trisurf(TR)

patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
hold on
% incenters = TR.incenter;
% faceNormals = TR.faceNormal;
% quiver3(incenters(:,1), incenters(:,2), incenters(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0)
view(3)

axis equal tight
