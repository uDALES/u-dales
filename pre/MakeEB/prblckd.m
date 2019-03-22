function [ prcntgblckd ] = prblckd(i,j,coa,ndima,sun,v1,cob,ndimb,F,centerweight,cornerweight,nblocks,nbw,bl)
%prcntgblckd     = percentage of view blocked
%i               = index of facet 1
%j               = index of facet 2      (only if sun==false)
%coa             = corners of facet 1
%ndima           = dimension of facet 1
%sun             = calculation if facet can see sun, not if two facets can
%                  see each other
%v1              = vector to the sun
%cob             = corners of facet 2    (only if sun==false)
%ndimb           = dimension of facet 2  (only if sun==false)
%F               = facets
%centerweight    = % of viewfield if centers see each other
%cornerweight    = % per corner if they see each other
%nblocks         = number of blocks
%nbw             = number of bounding walls
%bl              = list of blocks
testcrit=0; %=4;
prcntgblckd=0;
flag=0;
if sun  %calculation between facet 1 and the sun
    
    
    bi=F(i,3); %block index
    ci=F(i,4); %building index (-1 for roads, -99 for bounding wall)
    if ndima>testcrit  %also test corner, otherwise only test center
        
        for k=1:5
            facetpoint=coa(k,:);
            for n=1:nblocks %check if any block intersects
                if ((n==bi) && (ci>0))  %don't intersect with itself
                    continue
                end
                
                [flag,dint] = rbi(facetpoint, v1, bl(n,:));
                
                intersection = facetpoint + dint*v1;
                if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                    flag=0;
                end
                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                    if k==1
                        prcntgblckd=centerweight;
                    else
                        prcntgblckd=prcntgblckd+cornerweight;
                    end
                    break %get out of inner for loop
                end
            end
            
            if ~flag
                for m=1:nbw %check if any bounding wall intersects
                    
                    [flag,dint] = rbi(facetpoint, v1, bl(m+nblocks,:));
                    
                    intersection = facetpoint + dint*v1;
                    
                    if intersection(3)<facetpoint(3) %downstream direction of sun
                        flag=0;
                    end
                    if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                        if k==1 %it's the center of the facet
                            prcntgblckd=centerweight;
                        else
                            prcntgblckd=prcntgblckd+cornerweight;
                        end
                        break %get out of inner for loop
                    end
                end
            end
        end
    else
        % %disp('small area, checking only center')
        % %disp(['facet: ' num2str(i)])
        facetpoint=coa(1,:);  %center of facet
        for n=1:nblocks %check if any block intersects
            if ((n==bi) && (ci>0))  %don't intersect with itself
                continue
            end
            
            [flag,dint] = rbi(facetpoint, v1, bl(n,:));
            
            intersection = facetpoint + dint*v1;
            if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                % intersectionss(:)=0;
                flag=0;
            end
            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                prcntgblckd=1;
                break %get out of inner for loop
            end
        end
        
        if ~flag
            for k=1:nbw %check if any bounding wall intersects
                
                [flag,dint] = rbi(facetpoint, v1, bl(k+nblocks,:));
                
                intersection = facetpoint + dint*v1;
                
                if intersection(3)<facetpoint(3) %downstream direction of sun
                    flag=0;
                end
                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                    prcntgblckd=1;
                    break %get out of inner for loop
                end
            end
        end
        
    end
else %calculation between facet 1 and facet 2
    %disp('facet to facet')
    %disp(['% blocked so far ' num2str(prcntgblckd)])
    
    bi=F(i,3); %block index
    ci=F(i,4); %building index (-1 for roads, -99 for bounding wall)
    bi2=F(j,3); %block index
    fi=F(i,1); %facet index
    fi2=F(j,1); %facet index
    
    if (ndima+ndimb)>(2*testcrit)  %also test corner, otherwise only test center
        %disp('test corners')
        %disp(['% blocked so far ' num2str(prcntgblckd)])
        %find orientation to know which corners to compare
        ordera=[1 2 3 4 5];
        cas=fi*10+fi2;
        %disp('case')
        switch cas
            case 12, orderb=[1 3 4 5 2];
            case 13, orderb=[1 2 5 4 3];
            case 15, orderb=[1 3 2 5 4];
            case 21, orderb=[1 5 2 3 4];
            case 24, orderb=[1 5 4 3 2];
            case 31, orderb=[1 2 5 4 3];
            case 35, orderb=[1 5 4 3 2];
            case 42, orderb=[1 5 4 3 2];
            case 51, orderb=[1 3 2 5 4];
            case 53, orderb=[1 5 4 3 2];
            otherwise %14, 41, 23, 32, 25, 52, 34, 43, 45, 54
                orderb=[1 2 3 4 5];
        end
        
        for k=1:5
            test=1;
            oa=ordera(k);
            ob=orderb(k);
            facetpoint=coa(oa,:);
            facetpoint2=cob(ob,:);
            v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
            dint=norm(v);
            v1=v/dint;
            if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facets whos center is lower
                test=0;
            elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                test=0;
            elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                test=0;
            elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                test=0;
            elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                test=0;
            elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                test=0;
                % the following cases should be blocked by a block anyway, but
                % this is faster
            elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                test=0;
            elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                test=0;
            elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                test=0;
            elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                test=0;
            end
            if test
                %                 disp(['can potentially see each other' num2str(k)])
                %                 disp(['% blocked so far ' num2str(prcntgblckd)])
                for n=1:nblocks %check if any block intersects
                    if ((n==bi) && (ci>0))  %don't intersect with itself
                        continue
                    end
                    
                    [flag,dints] = rbi(facetpoint, v1, bl(n,:));
                    
                    if dints<=(0+eps)  %intersection is downstream
                        flag=0;
                    elseif dints>=(dint-eps) %block is behind target subfacet
                        flag=0;
                    end
                    if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                        if k==1
                            prcntgblckd=centerweight;
                        else
                            %disp(['blocked by ' num2str(n)])
                            prcntgblckd=prcntgblckd+cornerweight;
                        end
                        break %get out of inner for loop
                    end
                end
            else %facets can't possibly see each other
                %disp('facets cant possibly see each other')
                if k==1
                    prcntgblckd=centerweight;
                else
                    prcntgblckd=prcntgblckd+cornerweight;
                end
            end
        end
        
    else %small area, checking only center
        % %disp('small area, checkig only center')
        test=1;
        facetpoint=coa(1,:);
        facetpoint2=cob(1,:);
        v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
        dint=norm(v);
        v1=v/dint;
        
        %abs(v)
        %%disp(['v1: ' num2str(v1)])
        %test if facets can potentially see each other
        if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facets whos center is lower
            test=0;
        elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
            test=0;
        elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
            test=0;
        elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
            test=0;
        elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
            test=0;
        elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
            test=0;
            % the following cases should be blocked by a block anyway, but
            % this is faster
        elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
            test=0;
        elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
            test=0;
        elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
            test=0;
        elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
            test=0;
        end
        
        %%disp(['test: ' num2str(test)])
        
        if test %facets could see each other, check if blocked
            %calculate actual intersection distance
            for n=1:nblocks %check if any block intersects
                if ((n==bi) && (ci>0))  %don't intersect with itself
                    continue
                end
                
                [flag,dints] = rbi(facetpoint, v1, bl(n,:));
                
                if dints<=0  %intersection is downstream
                    flag=0;
                elseif dints>=dint %block is behind target subfacet
                    flag=0;
                end
                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                    %  %disp(['intersection with block: ' num2str(n)])
                    prcntgblckd=1;
                    break %get out of inner for loop
                end
            end
            
        else %facets can't possibly see each other
            prcntgblckd=1;
        end
    end
    
    
end