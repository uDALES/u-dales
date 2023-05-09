%% ULG that produces stl files, using road network method. 
% This code starts with lamp = 1 and subtracts until it hits its value
% First add green spaces us%% Define the input parameters

X = 96;                    % Domain length in computational points
Y  = 96;                   % Domain width in computational points
Z = 96;                    % Domain height in computatoinal points
zsize = 96;               % Height [m]
xsize = 96;               % Length [m]
ysize = 96;               % Width [m]
resz = zsize/Z;
resx = xsize/X;
resrat = resx/resz;
lp = 0.3;                   % Built fraction
lv = 0.1;                   % Vegetated fraction
lfr = 0.2;                  % Frontal aspect ratio we want 
lf = lfr*resrat;            % Accounting for the resolution 
minx = 4;                   % Narrowest a building can be in x
miny = 4;                   % Narrowest a building can be in y
minxg = 2;                  % Narrowest a green can be in x
minyg = 2;                  % Narrowest a green can be in y
minzr = 8;                  % Minimum building height [m]
maxzr = 30;                 % Max building height in [m]
minz  = round(minzr/resz);    % This stops buildings getting to small at fine resolutions
maxz = round(maxzr/resz);    % This stops buildings getting to tall at course resolutions
ngreen = 2;                 % Number of green spaces
max_chunk = 100;            % Limits how much the frontal surface area of building can increase at any one time
padding = 2;               % Minimum gap between buildings and edge
paddingg = 2;               % Minimum gap between green spaces and anything elseing the original approach
edgepaddingg = 4;           % Minimum gap between green and edge 
area   = X*Y;
builtarea = lp*area;
greenarea = lv*area;
front = lf*area;
expnr = '402'; 
tiled = false;
xtiles = 4;
ytiles = 2;


%% Make green patches
gnwidth = floor(sqrt(greenarea));
gnlength = floor(greenarea/gnwidth);
ptchs = [0,gnwidth,0,gnlength];
numptch = 1;
while numptch < ngreen
    ptchs = cutfromlist(ptchs, minxg,minyg);
    numptch = numptch + 1;
end    

ptchs_tall = ones(length(ptchs(:,1)),6);  %Prepare list in which green has height
ptchs_tall(:,1:4) = ptchs(:,1:4); 


%% Placing the green space randomly
topo = zeros(X, X);
first_topo = topo;       % Originally no blocks have been added
first_list = [];         % Originally no blocks have been added

% Add the green
[green_list, second_topo,green_outline] = placer(ptchs_tall,paddingg,edgepaddingg,X,Y,first_topo,first_list);

%% Add the buildings and reshape them to match lamp
new_blocks = blocks_based_on_green(X,Y,padding,green_list,ngreen);

changed_blocks = changing_lamp(new_blocks, lp,area,minx,miny,padding, green_list);
%changed_blocks = new_blocks

%% Give the buildings height

heights = heightadder2(changed_blocks,minz,maxz,front);

%% Check the frontal ratio is satisfied
total_front = sum(heights.*(changed_blocks(:,4)-changed_blocks(:,3)));
ratio = total_front/front;
final_heights = round(heights./ratio);
new_ratio = sum(final_heights.*(changed_blocks(:,4)-changed_blocks(:,3)))/front; % Rescale if need be


%% Add the building heights to the block list that describes them 
changed_blocks(:,6) = final_heights; % Add 1 to heights because the buildings start at z = 1 not z = 0
changed_blocks(:,5) = 0;

%% Bring in the green

final_list = [changed_blocks;green_list];

plotter(changed_blocks, green_list, X, Y, padding)
%% Make green areas the tallest to be consistent with uDALES preprocessing 
%maxh = max(final_list(:,6));
final_list(final_list(:,6) == 1,6) = 0; 
final_list(:,5) = 0; 
buildings = final_list;

%% Determine the actual lambda values we have achieved
lams = lambcalc(buildings,area,resrat);

%% Try transforming the coordinates to see if preprocessing params can match
% Here we have to add to the x and y lower coords and take away from the z
% upper coord such that the preprocessing recreates the geom we want.
% Note that these corrections are very different to those needed to convert
% from the geometries made in Python. I expect this is because Python
% indexes from 0 not 1...



if tiled == true
    repeated_buildings = [];
    repeated_green_outline = [];
    for xadd = 1:xtiles
        for yadd = 1:ytiles
            for b = 1:length(buildings(:,1))
                temp = [];
%                 new_xl = buildings(b,1)+(xadd-1)*X;
%                 new_xu = buildings(b,2)+(xadd-1)*X;
%                 new_yl = buildings(b,3)+(yadd-1)*Y;
%                 new_yu = buildings(b,4)+(yadd-1)*Y;
                new_xl = buildings(b,1)+xadd*X;
                new_xu = buildings(b,2)+xadd*X;
                new_yl = buildings(b,3)+yadd*Y;
                new_yu = buildings(b,4)+yadd*Y;
                temp = [new_xl, new_xu, new_yl,new_yu, buildings(b,5), buildings(b,6)];
                repeated_buildings = [repeated_buildings; temp];
            end 
            for g = 1:length(green_outline(:,1))
                temp = [];
%                 new_xl = green_outline(g,1)+(xadd-1)*X;
%                 new_xu = green_outline(g,2)+(xadd-1)*X;
%                 new_yl = green_outline(g,3)+(yadd-1)*Y;
%                 new_yu = green_outline(g,4)+(yadd-1)*Y;
                new_xl = green_outline(g,1)+xadd*X;
                new_xu = green_outline(g,2)+xadd*X;
                new_yl = green_outline(g,3)+yadd*Y;
                new_yu = green_outline(g,4)+yadd*Y;
                temp = [new_xl, new_xu, new_yl,new_yu];
                repeated_green_outline = [repeated_green_outline; temp];    
            end
        end
    end  
    buildings = repeated_buildings;
    green_outline = repeated_green_outline;
end 

function newblocks = cut(block, minx, miny) % Cuts a given block randomly i
% into 2 rectangles with minimum width and length
    disp('cutting')
    xl = 0;
    xu = block(2);
    yl = 0;
    yu = block(4);

    if xu < 2*minx & yu < 2*miny
        disp('too thin');
    elseif yu < 2*miny
        disp('too thin in y, cut along constant x');
        cutdir = 1;
    elseif xu < 2*minx    
        disp('too thin in x, force cut along constant y');
        cutdir = 2;
    else
        cutdir = randi(2);
    end 

    if cutdir == 1
        disp('cut along constant x, cutting in y direction')
        cutline = randi([minx,xu-minx]);
        newblock1 = [0,cutline,0,yu];
        newblock2 = [0,xu-cutline, 0,yu];
    elseif cutdir == 2
        disp('cut along constant y, cutting in x direction')
        cutline = randi([miny,yu-miny]);
        newblock1 = [0,xu,0,cutline];
        newblock2 = [0,xu, 0,yu- cutline];
    else 
        disp('too thing to cut')
    end 
    newblocks = [newblock1; newblock2];
end    

function newlist = cutfromlist(list, minx,miny) % Takes a list of blocks,
% picks one and cuts it, picks new block if current is too small to cut 
    n = size(list,1);
    cutind = randi([1,n]);
    block2cut = list(cutind,:);
    xu = block2cut(2);
    yu = block2cut(4);
    while xu < 2*minx && yu<2*miny
        disp('block too small to cut in x or y')
        cutind = randi([1,n]);
        block2cut = list(cutind,:);
        xu = block2cut(2);
        yu = block2cut(4);
    end    
    list(cutind,:) = [];
    newblocks = cut(block2cut, minx,miny);
    newlist = [list;newblocks];
end    

function [new_list, new_topo, outline] = placer(blocks,padding,edgepadding,X,Y,old_topo,old_list)
% Takes a set of blocks, and places them randomly on a pre-existing
% arrangement, each block being no closer to another than a padding
% distance. It then returns a list consisting of the original blocks and
% the new ones that have been added (the list is in [xl,xu,yl,yu,zl,zu;
% ...] form. It also returns a new topo with all the blocks in it as well
% as an outline, giving the locations of the edges of the new blocks that
% were added.
    blocks(:,1:4) = blocks(:,1:4)+1; %To avoid indexing from zero issues
    empty = zeros(Y,X);
    outline = [];
    i = 1;
    topo = old_topo; 
    while i <= length(blocks(:,1))
        attempts = 0;
        blk = blocks(i,:);
        xlmax = X - (blk(2)+edgepadding);
        ylmax = Y - (blk(4)+edgepadding);
        fit = false;
        while fit == false
            attempts = attempts + 1;
            if attempts > 10000
                i = 1;
                topo = old_topo;
                outline = [];
                break
            end
            trial = empty;
            xadd = randi([edgepadding,xlmax]);
            yadd = randi([edgepadding,ylmax]);
            shft(1) = blk(1)+xadd;
            shft(2) = blk(2)+xadd;
            shft(3) = blk(3)+yadd;
            shft(4) = blk(4)+yadd;
            trial(shft(3):shft(4),shft(1):shft(2)) = blk(5);
            padded = empty;
            padded(shft(3)-padding:shft(4)+padding,shft(1)-padding:shft(2)+padding) = blk(5);
            if (padded.*topo == 0)
                topo = topo + trial;
                old_list = [old_list;shft(1),shft(2),shft(3), shft(4), 1, blk(5)];
                fit = true;
                outline = [outline; shft(1),shft(2),shft(3), shft(4)];
                i = i+1;
            end
        end
    end 
    new_list = old_list;
    new_topo = topo;
end

function building_blocks = blocks_based_on_green(X,Y,padding,green_list, ngreen)
    % Cutting out space based on location of green corners
    xs = sort(unique([padding, X-padding,reshape(green_list(:,1:2),1,numel(green_list(:,1:2)))]));
    ys = sort(unique([padding,Y-padding, reshape(green_list(:,3:4),1,numel(green_list(:,3:4)))]));
    
    new_blocks = [];
    nblocks_removed = 0;
    
    for i =1:length(xs)-1
         for j = 1:length(ys)-1
             new_blocks = [new_blocks; xs(i)+padding, xs(i+1)-padding, ys(j)+padding, ys(j+1)-padding,1,1];
         end
    end
    for i = length(new_blocks(:,1)):-1:1
        if or(new_blocks(i,1) >= new_blocks(i,2), new_blocks(i,4) <= new_blocks(i,3))
            new_blocks(i,:) = [];   
        end
    end
    original_newblocks = new_blocks;
    removed_blocks = [];
    %% We need to remove blocks that correspond to the green blocks that are in the list, new_blocks
    for g = 1:length(green_list(:,1))
        gxl = green_list(g,1);
        gxu = green_list(g,2);
        gyl = green_list(g,3);
        gyu = green_list(g,4);
        polyg = polyshape([gxl gxl gxu gxu],[gyl gyu gyu gyl]);
        for i = length(new_blocks(:,1)):-1:1
            bxl = new_blocks(i,1);
            bxu = new_blocks(i,2);
            byl = new_blocks(i,3);
            byu = new_blocks(i,4);
            polyb = polyshape([bxl bxl bxu bxu],[byl byu byu byl]);
            polytest = intersect(polyb, polyg);
            inter = area(polytest);
            if inter ~= 0 
                disp(inter)
                disp('Green to be removed.')
                removed_blocks = [removed_blocks; new_blocks(i,:)];
                new_blocks(i,:) = [];
                nblocks_removed = nblocks_removed + 1;
                if nblocks_removed > ngreen
                    disp('wtf')
                end 
            end     
        end
    end
    building_blocks = new_blocks;
end

function new_heights = heightadder2(blocks, minz,maxz, front)
% Here we pick a random building then add a random height to it. 
    side_lens = blocks(:,4)-blocks(:,3);
    %tot_len = sum(side_lens);
    min_areas = side_lens*minz;                     % Frontal areas based on minimum heights
    remaining = front - sum(min_areas);             % How much area we need to add until lf is satisfied
    heights = minz*ones(length(blocks(:,1)),1);     % Each block is initially its minimum height
    blockind_list = 1:length(blocks(:,4));
    while remaining > 0
        if length(blockind_list) == 0
            disp('No more blocks we can add to.')
            break
        end     
        if remaining < min((blocks(:,4)-blocks(:,3)))
            disp('Insufficient area remaining to add 1m to smallest building.')
            break
        end    
        blockind_ind = randi(length(blockind_list));
        blockind = blockind_list(blockind_ind);

        max_add = maxz-heights(blockind);
        if max_add <= 0 
            disp('block max height')
            blockind_list(blockind_list == blockind) = [];
            continue
        end
        blockwidth = (blocks(blockind,4)- blocks(blockind,3));
        max_add = min(max_add,floor(remaining/blockwidth));
        if max_add <=0      
            disp('Not enough remaining to add to this block')
            blockind_list(blockind_list == blockind) = [];
            continue
        end    
        height_add = randi(max_add);
        chunk = height_add*blockwidth;
        heights(blockind) = heights(blockind)+height_add;
        remaining = remaining - chunk;
    end
    new_heights = heights;
end

function new_blocks = changing_lamp(blocks,lp,Area, minx,miny,padding, green_list)
        areas = (blocks(:,2)- blocks(:,1)).*(blocks(:,4)- blocks(:,3));
        
%         %% Remove blocks that are too small
%         for i = length(blocks(:,1)):-1:1
%             if or(blocks(i,2)- blocks(i,1)<minx, blocks(i,4)- blocks(i,3)<miny) 
%                 Area1 = areas(i);
%                 disp('Thin block removed.')
%                 blocks(i,:) = [];
%                 areas(i) = [];
%                 remaining = remaining - Area1;
%             end
%         end
        lamp = sum(areas)/Area;
        if lamp > lp
            disp('More must be removed.')
            remaining = sum(areas)-round(lp*Area);
            while remaining>0.01*lp*Area
                %length(blocks(:,1))
                for i = length(blocks(:,1)):-1:1
                    if or(blocks(i,2)- blocks(i,1)<minx, blocks(i,4)- blocks(i,3)<miny) 
                        Area1 = areas(i);
                        if Area1 <= remaining
                            disp('Thin block removed.')
                            blocks(i,:) = [];
                            remaining = remaining - Area1;
                        end
                    end
                end
                for i = length(blocks(:,1)):-1:1 
                    if and(blocks(i,2)- blocks(i,1)<minx,blocks(i,4)- blocks(i,3)>=miny)
                        disp('Thin block too thin in x, cut in y ')
                        len = blocks(i,2)- blocks(i,1);
                        best_cut = round(remaining/len);
                        new_yu = max(blocks(i,3)+miny,blocks(i,4) - best_cut);
                        cut = blocks(i,4) - new_yu;
                        blocks(i,4) = new_yu;
                        remaining = remaining - (blocks(i,2)- blocks(i,1))*cut;
                    elseif and(blocks(i,2)- blocks(i,1)>=minx,blocks(i,4)- blocks(i,3)<miny)
                        disp('Thin block too thin in y, cut in x ')
                        width = blocks(i,4)- blocks(i,3);
                        best_cut = round(remaining/width);
                        new_xu = max(blocks(i,1)+minx,blocks(i,2) - best_cut);
                        cut = blocks(i,2) - new_xu;
                        blocks(i,2) = new_xu;
                        remaining = remaining - cut*(blocks(i,4)-blocks(i,3));
                    end
                end
                while remaining > 0.01*lp*Area
                    if remaining < minx
                        break
                    end    
                    for i = 1:length(blocks(:,1))
                        if and(blocks(i,2)- blocks(i,1)>minx,blocks(i,4)- blocks(i,3)>miny)
                            disp('Thick enough in both directions')
                            if blocks(i,2)- blocks(i,1) < blocks(i,4)- blocks(i,3)
                                disp('Block longer in y so reduce yu')
                                len = blocks(i,2)- blocks(i,1);
                                best_cut = round(remaining/len);
                                %new_yu = max(blocks(i,3)+miny,blocks(i,4) - best_cut);
                                %new_yu = randi([min(blocks(i,3)+miny,blocks(i,4) - best_cut), max(blocks(i,3)+miny,blocks(i,4) - best_cut)]);
                                new_yu = randi([blocks(i,3)+miny,blocks(i,4)]);
                                cut = blocks(i,4) - new_yu;
                                if remaining - (blocks(i,2)- blocks(i,1))*cut > - 0.01*lp*Area
                                    blocks(i,4) = new_yu;
                                    remaining = remaining - (blocks(i,2)- blocks(i,1))*cut;
                                end 
                            elseif blocks(i,2)- blocks(i,1) > blocks(i,4)- blocks(i,3)
                                disp('Block longer in x so reduce xu')
                                width = blocks(i,4)- blocks(i,3);
                                best_cut = round(remaining/width);
                                %new_xu = max(blocks(i,1)+minx,blocks(i,2) - best_cut);
                                new_xu = randi([blocks(i,1)+minx, blocks(i,2)]);
                                cut = blocks(i,2) - new_xu;
                                if remaining - cut*(blocks(i,4)-blocks(i,3)) > - 0.01*lp*Area
                                    blocks(i,2) = new_xu;
                                    remaining = remaining - cut*(blocks(i,4)-blocks(i,3))
                                end 
                            else     
                                disp('Block equal in x and y so reduce xu')
                                width = blocks(i,4)- blocks(i,3);
                                best_cut = round(remaining/width);
                                new_xu = randi([blocks(i,1)+minx,blocks(i,2)]);
                                cut = blocks(i,2) - new_xu;
                                if remaining - cut*(blocks(i,4)-blocks(i,3)) > - 0.01*lp*Area
                                    blocks(i,2) = new_xu;
                                    remaining = remaining - cut*(blocks(i,4)-blocks(i,3));
                                end 
                            end
                        end
                    end     
                                


                end
            end   

        elseif lamp == lp
            disp('Lam is correct.')
        else 
            disp('More must be added.')
            remaining = round(lp*Area)- sum(areas);
            %num_blocks = length(blocks(:,1))
            xmerge = true;
            ymerge = false; 
            attempts = 1000;
            while attempts > 0
                if remaining > 0.01*lp*Area
                        if xmerge 
                        [blocks, remaining] = mergerx(blocks, remaining, green_list, Area, lp);
                        attempts = attempts -1;
                        xmerge = false;
                        ymerge = true;
                        disp('xmerg');
                        elseif ymerge
                            [blocks, remaining] = mergery(blocks, remaining, green_list, Area, lp);
                            attempts = attempts -1;
                            xmerge = true;
                            ymerge = false; 
                            disp('ymerg');
                        end    
                %num_blocks = length(blocks(:,1));
                else 
                    attempts = 0;
                end
            end
            new_blocks = blocks;
        end    
        new_blocks = blocks;
end        

function lams = lambcalc(blocks,Area, resrat)
    buildings = blocks(blocks(:,6) ~= max(blocks(:,6)),:); 
    builtareas = (buildings(:,2)-buildings(:,1)).*(buildings(:,4)-buildings(:,3));
    builtarea = sum(builtareas);
    lamp = builtarea/Area;

    frontareas = (buildings(:,4)-buildings(:,3)).*(buildings(:,6)-buildings(:,5));
    frontarea = sum(frontareas);
    lamf = frontarea/Area;

    greens = blocks(blocks(:,6) == max(blocks(:,6)),:);
    greenareas = (greens(:,2)-greens(:,1)).*(greens(:,4)-greens(:,3));
    greenarea = sum(greenareas);
    lamv = greenarea/Area;

    lams = [lamp,lamv,lamf/resrat];
end

function [new_blocks, remains] = mergerx(blocks, remaining, green_list, Area, lp)
    i = randi(length(blocks(:,1)));
    block = blocks(i,:);
    otherblocks = blocks;
    otherblocks(i,:) = [];
    xl = block(1);
    xu = block(2);
    yl = block(3);
    yu = block(4);
    zl = block(5);
    zu = block(6);
    merg_pos = [];
    for b = 1:length(otherblocks(:,1))
        checkblock = otherblocks(b,:);
        cxl = checkblock(1);
        cxu = checkblock(2);
        cyl = checkblock(3);
        cyu = checkblock(4);
        if and(yl == cyl, yu == cyu)
            %disp('merge possible')
            merg_pos = [merg_pos; checkblock]; % Populate array of blocks that are at the right y location
        end
    end
    if length(merg_pos) == 0
        disp('No merge possible')
        new_blocks = blocks;
        remains = remaining;
        return
    end
    min_gaps = min(abs(xl-merg_pos(:,2)),abs(xu-merg_pos(:,1)) ); % Find the smallest gap between this block and the others in line with it
    min_gap = min(min_gaps); % Identify which block is the closest 
    merg_candidate_index = find(min_gaps == min_gap);
    if (length(merg_candidate_index) == [])
        disp('no remaining merg')
    elseif (length(merg_candidate_index)>1)
        merg_candidate_index = merg_candidate_index(1);
    end  
    merg_candidate = merg_pos(merg_candidate_index,:);
    
    new_block = [min(xl,merg_candidate(1)),max(xu,merg_candidate(2)),yl,yu,1,1];
    nxl = new_block(1);
    nxu = new_block(2);
    nyl = new_block(3);
    nyu = new_block(4);
    npoly = polyshape([nxl,nxl,nxu,nxu],[nyl,nyu,nyu,nyl]);
    gn_overlap = 0;
     
    for g = 1:length(green_list(:,1))
        gxl = green_list(g,1);
        gxu = green_list(g,2);
        gyl = green_list(g,3);
        gyu = green_list(g,4);
        gpoly = polyshape([gxl,gxl,gxu,gxu],[gyl,gyu,gyu,gyl]);
        gn_overlap = gn_overlap + area(intersect(npoly,gpoly));
    end
    if gn_overlap == 0
        
        %disp('Merg possible')
        gry_overlap = 0;
        for gy = 1:length(otherblocks(:,1))
            gyxl = otherblocks(gy,1);
            gyxu = otherblocks(gy,2);
            gyyl = otherblocks(gy,3);
            gyyu = otherblocks(gy,4);
            gypoly = polyshape([gyxl,gyxl,gyxu,gyxu],[gyyl,gyyu,gyyu,gyyl]);
            overlap = area(intersect(npoly,gypoly));
            if overlap > 0 
                gry_overlap = gry_overlap + 1;
            end 
        end  
        if gry_overlap > 1
            disp('overlap')
            new_blocks = blocks;
        else 
            gain = (yu-yl)*min_gap;
            if (remaining-gain)/(lp*Area) < -0.01
                disp('section is too big to add')
            else
                disp('Suitable to remove')
                blocks(i,:) = [];
                [~, block2go] = ismember(merg_candidate, blocks,'rows');
                blocks(block2go,:) = [];         
                disp('Original blocks removed.')
                blocks = [blocks ; new_block];
                disp('New block added.')
                remaining = remaining - gain;
                gain
                disp('MERGE!!!!!!!');
            end 
        end 
        
    else 
        disp('No merg possible')
        new_blocks = blocks;
    end
    remains = remaining;
    new_blocks = blocks;
end 

function [new_blocks, remains] = mergery(blocks, remaining, green_list, Area, lp)
    i = randi(length(blocks(:,1)));
    block = blocks(i,:);
    otherblocks = blocks;
    otherblocks(i,:) = [];
    xl = block(1);
    xu = block(2);
    yl = block(3);
    yu = block(4);
    zl = block(5);
    zu = block(6);
    merg_pos = [];
    for b = 1:length(otherblocks(:,1))
        checkblock = otherblocks(b,:);
        cxl = checkblock(1);
        cxu = checkblock(2);
        cyl = checkblock(3);
        cyu = checkblock(4);
        if and(xl == cxl, xu == cxu)
            %disp('merge possible')
            merg_pos = [merg_pos; checkblock]; % Populate array of blocks that are at the right y location
        end
    end
    if length(merg_pos) == 0
        disp('No merge possible')
        new_blocks = blocks;
        remains = remaining;
        return
    end
    min_gaps = min(abs(yl-merg_pos(:,4)),abs(yu-merg_pos(:,3)) ); % Find the smallest gap between this block and the others in line with it
    min_gap = min(min_gaps); % Identify which block is the closest 
    merg_candidate_index = find(min_gaps == min_gap);
    if (length(merg_candidate_index) == [])
        disp('no remaining merg')
    elseif (length(merg_candidate_index)>1)
        merg_candidate_index = merg_candidate_index(1);
    end  
    merg_candidate = merg_pos(merg_candidate_index,:);
    
    new_block = [xl,xu,min(yl,merg_candidate(3)),max(yu,merg_candidate(4)),1,1];
    nxl = new_block(1);
    nxu = new_block(2);
    nyl = new_block(3);
    nyu = new_block(4);
    npoly = polyshape([nxl,nxl,nxu,nxu],[nyl,nyu,nyu,nyl]);
    gn_overlap = 0;
     
    for g = 1:length(green_list(:,1))
        gxl = green_list(g,1);
        gxu = green_list(g,2);
        gyl = green_list(g,3);
        gyu = green_list(g,4);
        gpoly = polyshape([gxl,gxl,gxu,gxu],[gyl,gyu,gyu,gyl]);
        gn_overlap = gn_overlap + area(intersect(npoly,gpoly));
    end
    if gn_overlap == 0
        
        %disp('Merg possible')
        gry_overlap = 0;
        for gy = 1:length(otherblocks(:,1))
            gyxl = otherblocks(gy,1);
            gyxu = otherblocks(gy,2);
            gyyl = otherblocks(gy,3);
            gyyu = otherblocks(gy,4);
            gypoly = polyshape([gyxl,gyxl,gyxu,gyxu],[gyyl,gyyu,gyyu,gyyl]);
            overlap = area(intersect(npoly,gypoly));
            if overlap > 0 
                gry_overlap = gry_overlap + 1
            end 
        end  
        if gry_overlap > 1
            disp('overlap')
            new_blocks = blocks;
        else 
            gain = (xu-xl)*min_gap
            if (remaining-gain)/(lp*Area) < -0.01
                disp('section is too big to add')
            else
                disp('Suitable to remove')
                blocks(i,:) = [];
                [~, block2go] = ismember(merg_candidate, blocks,'rows');
                blocks(block2go,:) = [];         
                disp('Original blocks removed.')
                blocks = [blocks ; new_block];
                disp('New block added.')
                remaining = remaining - gain;
                disp('MERGE!!!!!!!');
            end 
        end 
        
    else 
        disp('No merg possible')
        new_blocks = blocks;
    end
    remains = remaining;
    new_blocks = blocks;    
end 

function plotter(bs,gs,xlim,ylim,padding)
xl = 0; 
xu = xlim;
yl = 0;
yu = ylim;
x = [xl, xu, xu, xl];
y = [yl, yl, yu, yu];
z = [0, 0, 0, 0];
fh = figure;
hold on
light();
clr = [0.0,0.0,0.0];
patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7); 
view(3)
axis equal;
axis off;
grid off;
ax = gca;
ax.FontSize = 18;
xlabel('$x(\mathrm{m})$','interpreter','latex', 'FontSize',22);
ylabel('$y(\mathrm{m})$','interpreter','latex', 'FontSize',22);
zlabel('$z(\mathrm{m})$','interpreter','latex', 'FontSize',22);
clr = [0.0,0.75,0.0];
z = [0.01,0.01,0.01,0.01];
for g = 1:length(gs(:,1))
    xl = gs(g,1)+padding;
    xu = gs(g,2)-padding;
    yl = gs(g,3)+padding;
    yu = gs(g,4)-padding;
    x = [xl, xu, xu, xl];
    y = [yl, yl, yu, yu];
    patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
end 
clr = [0.75,0.75,0.75];
hold on
    for b = 1:length(bs(:,1))
        xl = bs(b,1);
        xu = bs(b,2);
        yl = bs(b,3);
        yu = bs(b,4);
        zl = 0;
        zu = bs(b,6);
        %% top
        x = [xl, xu, xu, xl];
        y = [yl, yl, yu, yu];
        z = [zu, zu, zu, zu];
        patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
        %% west
        x = [xl, xl, xl, xl];
        y = [yl, yu, yu, yl];
        z = [zu, zu, zl, zl];
        patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
        %% east
        x = [xu, xu, xu, xu];
        y = [yl, yu, yu, yl];
        z = [zu, zu, zl, zl];
        patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
        %% north
        x = [xl, xu, xu, xl];
        y = [yu, yu, yu, yu];
        z = [zu, zu, zl, zl];
        patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
        %% south 
        x = [xl, xu, xu, xl];
        y = [yl, yl, yl, yl];
        z = [zu, zu, zl, zl];
        patch('XData',x,'YData',y,'ZData',z,'FaceColor',clr,'edgecolor','none','FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
    end
    %pause(0.1) 
    %set(gca,'View',[180 90]);
end    