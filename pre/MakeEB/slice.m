function [coaa, cobb, pf1sf2u] = slice(fi,fi2,coa,cob,cornerweight,centerweight)
%in case one facet cuts through
%the plane of the other facet:
%                    |          .-----------.
%   ---------------  |    .or.  |           |
%                    |          | . . . . . |.------.
%                    |          |___________||      |
%                                            |______|
%
% if this is the case slice facet up (max 2 slices) -> coaa/cobb
% count the number of slices -> m
% increase the percentage they see of each other, since the two corners were blocked but now
% they are not included in the vf-calculation anymore -> pf1sf2u


%1=top,2=west face,3=east face,4=north face, 5=south face
pf1sf2u=0;
coaa=coa(6:9,:);
cobb=cob(6:9,:);
nottested=1; %necessary to check for i and j (alternatively the if's in each case could be combined)
if ((fi*fi2==6) || (fi*fi2==20)) %opposite each other, do nothing;
    return
elseif (fi==1) %fi2=2,3,4,5
    
    if (cobb(1,3)<coaa(1,3)) %j starts lower than height of i
        cobb(1,3)=coaa(1,3); %slice j on same height as i
        cobb(4,3)=coaa(1,3); %slice j on same height as i
        pf1sf2u=2*cornerweight; %add 2*cornerweight since lower corner now is visible
        nottested=0;
        if  cob(1,3)<=coaa(1,3) %if center of j is lower than height also add centerweight
            pf1sf2u=2*cornerweight+centerweight;
        end
    end
    
elseif (fi==2) %fi2=1,4,5
    
    
    if (((fi2==1) && (coaa(1,1)<cobb(4,1))) || (((fi2==4) || (fi2==5))&& (coaa(1,1)<cobb(4,1))))%my x<xmax
        cobb(3,1)=coaa(1,1);
        cobb(4,1)=coaa(1,1);
        pf1sf2u=2*cornerweight;
        nottested=0;
        if coaa(1,1)<=cob(1,1) %if j center is east of i wall
            pf1sf2u=2*cornerweight+centerweight  ;
        end
        % elseif (((fi2==4) || (fi2==5))&& (coa(1,1)<cob(4,1)))
        %     cobb=cob(2:5,:);
        %     cobb(3,1)=coaa(1,1);
        %     cobb(4,1)=coaa(1,1);
        %     pf1sf2u=2*cornerweight;
    end
elseif (fi==3) %fi2=1,4,5
    
    
    if (((fi2==1) && (coaa(1,1)>cobb(1,1))) || (((fi2==4) || (fi2==5)) && (coaa(1,1)>cobb(1,1)))) %my x>xmin
        cobb(1,1)=coaa(1,1);
        cobb(2,1)=coaa(1,1);
        pf1sf2u=2*cornerweight;
        nottested=0;
        if coaa(1,1)>=cob(1,1) %if j center is west of i wall
            pf1sf2u=2*cornerweight+centerweight  ;
        end
        % elseif (((fi2==4) || (fi2==5)) && (coa(1,1)>cob(1,1)))
        %     cobb=cob(2:5,:);
        %     cobb(1,1)=coaa(1,1);
        %     cobb(2,1)=coaa(1,1);
        %     pf1sf2u=2*cornerweight;
    end
    
elseif (fi==4) %fi2=1,2,3
    
    
    if ((fi2==1) && (coaa(1,2)>cobb(1,2))) %my y>ymin
        
        cobb(1,2)=coaa(1,2);
        cobb(4,2)=coaa(1,2);
        pf1sf2u=2*cornerweight;
        nottested=0;
        if coaa(1,2)>=cob(1,2) %if j center is south of i wall
            pf1sf2u=2*cornerweight+centerweight  ;
        end
    elseif (((fi2==2) || (fi2==3)) && (coaa(1,2)>cobb(1,2)))
        
        cobb(1,2)=coaa(1,2);
        cobb(2,2)=coaa(1,2);
        pf1sf2u=2*cornerweight;
        nottested=0;
        if coaa(1,2)>=cob(1,2) %if j center is south of i wall
            pf1sf2u=2*cornerweight+centerweight  ;
        end
    end
    
elseif (fi==5) %fi2=1,2,3
    %disp('fi==5')   
    %disp(['nottested: ' num2str(nottested)])

    if ((fi2==1) && (coaa(1,2)<cobb(2,2))) %my y<ymax
        cobb(2,2)=coaa(1,2);
        cobb(3,2)=coaa(1,2);
        pf1sf2u=2*cornerweight;
        nottested=0;
        if coaa(1,2)<=cob(1,2) %if j center is north of i wall
            pf1sf2u=2*cornerweight+centerweight  ;
        end
    elseif (((fi2==2) || (fi2==3)) && (coaa(1,2)<cobb(3,2))) %ILS13 11.11.17, was cobb(2,2)
        cobb(3,2)=coaa(1,2);
        cobb(4,2)=coaa(1,2);
        pf1sf2u=2*cornerweight;
        nottested=0;
        if coaa(1,2)<=cob(1,2) %if j center is north of i wall
            pf1sf2u=2*cornerweight+centerweight  ;
        end
    end
end
if nottested
    if (fi2==1)  %fi=2,3,4,5
        if (coaa(1,3)<cobb(1,3)) %i starts lower than height of j
            coaa(1,3)=cobb(1,3); %slice i on same height as j
            coaa(4,3)=cobb(1,3); %slice i on same height as j
            pf1sf2u=2*cornerweight; %add 2*cornerweight since lower corner now is visible
            if  coa(1,3)<=cobb(1,3) %if center of i is lower than height also add centerweight
                pf1sf2u=2*cornerweight+centerweight;
            end
        end
        %end
        
    elseif (fi2==2) %fi=1,4,5
        if (((fi==1) && (cobb(1,1)<coaa(4,1))) || (((fi==4) || (fi==5))&& (cobb(1,1)<coaa(4,1))))%my x<xmax
            coaa(3,1)=cobb(1,1);
            coaa(4,1)=cobb(1,1);
            pf1sf2u=2*cornerweight;
            
            if cobb(1,1)<=coa(1,1) %if i center is east of j wall
                pf1sf2u=2*cornerweight+centerweight  ;
            end
            
        end
        
    elseif (fi2==3) %fi=1,4,5
        if (((fi==1) && (cobb(1,1)>coaa(1,1))) || (((fi==4) || (fi==5)) && (cobb(1,1)>coaa(1,1)))) %my x>xmin
            coaa(1,1)=cobb(1,1);
            coaa(2,1)=cobb(1,1);
            pf1sf2u=2*cornerweight;
            if cobb(1,1)>=coa(1,1) %if i center is west of j wall
                pf1sf2u=2*cornerweight+centerweight  ;
            end
            
        end
        
    elseif (fi2==4) %fi=1,2,3
        if ((fi==1) && (cobb(1,2)>coaa(1,2))) %my y>ymin
            coaa(1,2)=cobb(1,2);
            coaa(4,2)=cobb(1,2);
            pf1sf2u=2*cornerweight;
            if cobb(1,2)>=coa(1,2) %if i center is south of j wall
                pf1sf2u=2*cornerweight+centerweight  ;
            end
        elseif (((fi==2) || (fi==3)) && (cobb(1,2)>coaa(1,2)))
            coaa(1,2)=cobb(1,2);
            coaa(2,2)=cobb(1,2);
            pf1sf2u=2*cornerweight;
            if cobb(1,2)>=coa(1,2) %if i center is south of j wall
                pf1sf2u=2*cornerweight+centerweight  ;
            end
        end
    elseif (fi2==5)  %fi=1,2,3
        if ((fi==1) && (cobb(1,2)<coaa(2,2))) %my y<ymax
            coaa(2,2)=cobb(1,2);
            coaa(3,2)=cobb(1,2);
            pf1sf2u=2*cornerweight;
            if cobb(1,2)<=coa(1,2) %if j center is north of i wall
                pf1sf2u=2*cornerweight+centerweight  ;
            end
        elseif (((fi==2) || (fi==3)) && (cobb(1,2)<coaa(2,2)))
            coaa(3,2)=cobb(1,2);
            coaa(4,2)=cobb(1,2);
            pf1sf2u=2*cornerweight;
            if cob(1,2)<=cob(1,2) %if j center is north of i wall
                pf1sf2u=2*cornerweight+centerweight  ;
            end
        end
    end
end
