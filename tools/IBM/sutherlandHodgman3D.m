%The inputs are a table of x-y pairs for the verticies of the subject
%polygon and boundary polygon. (x values in column 1 and y values in column
%2) The output is a table of x-y pairs for the clipped version of the 
%subject polygon.

function clippedPolygon = sutherlandHodgman3D(subjectPolygon,clipPlanes)

%% Helper Functions

    %computeIntersection() assumes the two lines intersect
    function intersection = computeIntersection(plane,line)
        % Line described by p1 + u(p2 - p1) = 0
        % Plane described by n . (p - p3) = 0
        p1 = line(1,:);
        p2 = line(2,:);
        p3 = plane(1:3) * plane(4);
        u = (plane(1:3) * (p3 - p1)') / (plane(1:3) * (p2 - p1)');
        intersection = p1 + u * (p2 - p1);

    end %computeIntersection

    function in = inside(point,plane)
        
        if (dot(point, plane(1:3)) - plane(4) <= 0)
            in = true;
        else
            in = false;
        end
        
    end %inside

%% Sutherland-Hodgman Algorithm

    clippedPolygon = subjectPolygon;
    %numVerticies = size(clipPolygon,1);
    numPlanes = size(clipPlanes);
    %clipVertexPrevious = clipPolygon(end,:);
    
    %for clipVertex = (1:numVerticies)
    for p = 1:numPlanes
        %clipBoundary = [clipPolygon(clipVertex,:) ; clipVertexPrevious]; % Need to be a plane
        clipPlane = clipPlanes(p,:);
        inputList = clippedPolygon;
        
        clippedPolygon = [];
        if ~isempty(inputList)
            previousVertex = inputList(end,:);
        end

        for subjectVertex = (1:size(inputList,1))
            if ( inside(inputList(subjectVertex,:),clipPlane) )
                if( not(inside(previousVertex,clipPlane)) )  
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    intersection = computeIntersection(clipPlane,subjectLineSegment);
                    %clippedPolygon(end+1,:) = computeIntersection(clipPlane,subjectLineSegment);
                    if isempty(clippedPolygon)
                        clippedPolygon(end+1,:) = intersection;
                    else
                        if ~any(ismember(clippedPolygon, intersection, 'rows'))
                            clippedPolygon(end+1,:) = intersection;
                        end
                    end
                end

                if isempty(clippedPolygon)
                    clippedPolygon(end+1,:) = inputList(subjectVertex,:);
                else
                    if ~any(ismember(clippedPolygon, inputList(subjectVertex,:), 'rows'))
                        clippedPolygon(end+1,:) = inputList(subjectVertex,:);
                    end
                end
                
            elseif( inside(previousVertex,clipPlane) )
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    intersection = computeIntersection(clipPlane,subjectLineSegment);
                    %clippedPolygon(end+1,:) = computeIntersection(clipPlane,subjectLineSegment);

                    if isempty(clippedPolygon)
                        clippedPolygon(end+1,:) = intersection;
                    else
                        if ~any(ismember(clippedPolygon, intersection, 'rows'))
                            clippedPolygon(end+1,:) = intersection;
                        end
                    end
            end
            
            previousVertex = inputList(subjectVertex,:);
            %clipVertexPrevious = clipPolygon(clipVertex,:);
            
        end %for subject verticies                
    end %for boundary verticies
end %sutherlandHodgman