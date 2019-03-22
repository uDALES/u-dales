function [ ndim, area, co] = detsub(i,F,BB,G,W,cornm,xb,yb,zb,delta)
% returns ndim=(dim1*dim2), area, center & corner of every facet (could be preprocessed and
% saved)

bi=F(i,3); %block index
fi=F(i,1); %facet index
ci=F(i,4); %building index (-1 for roads, -99 for bounding wall)

il=-1; iu=-1;jl=-1;ju=-1;kl=-1;ku=-1;it=-1;jt=-1;kt=-1;
co=-1;
ndim=-1;
area=-1;




if (ci<=-100) %it is a bounding wall
    kl=W(bi,5)+1;
    ku=W(bi,6)+1;
    if (fi==2)  %east, facing west
        jl=W(bi,3);     %lower y index of floor facet 1
        ju=W(bi,4);     %upper y index of floor facet 1
        it=W(bi,1);
        
        x=xb(it+1);
        
        ndim=((ju-jl)+1)*((ku-kl)+1);
        
    elseif (fi==3)  %west, facing east
        jl=W(bi,3);     %lower y index of floor facet 1
        ju=W(bi,4);     %upper y index of floor facet 1
        it=W(bi,2);
        
        x=xb(it);
        
        ndim=((ju-jl)+1)*((ku-kl)+1);
    elseif (fi==4)  %south, facing north
        il=W(bi,1);
        iu=W(bi,2);
        jt=W(bi,4);
        
        y=yb(jt);
        ndim=((iu-il)+1)*((ku-kl)+1);
        
    else %if (fi==5)  %north, facing south
        jt=W(bi,3);
        il=W(bi,1);
        iu=W(bi,2);
        
        y=yb(jt+1);
        
        ndim=((iu-il)+1)*((ku-kl)+1);
    end
elseif (ci>=0) %deal with floors separately
    if (fi==1)  %top
        il=BB(bi,1);     %lower x index of facet 1
        iu=BB(bi,2);     %upper x index of facet 1
        jl=BB(bi,3);     %lower y index of facet 1
        ju=BB(bi,4);     %upper y index of facet 1
        kt=BB(bi,6)+1;
        z=zb(kt+1);
        
        ndim=((iu-il)+1)*((ju-jl)+1);
        
    elseif ((fi==2) || (fi==3))  %west / east
        jl=BB(bi,3);     %lower y index of facet 1
        ju=BB(bi,4);     %upper y index of facet 1
        kl=BB(bi,5)+1;     %lower z index of facet 1
        ku=BB(bi,6)+1;     %upper z index of facet 1
        
        if fi==2  %west, facing west
            it=BB(bi,1);
            x=xb(it);
        else %east, facing east
            it=BB(bi,2);
            x=xb(it+1);
        end
        
        ndim=((ju-jl)+1)*((ku-kl)+1);
        
        
    elseif ((fi==4) || (fi==5)) %north / south
        il=BB(bi,1) ;    %lower x index of facet 1
        iu=BB(bi,2) ;    %upper x index of facet 1
        kl=BB(bi,5)+1 ;    %lower y index of facet 1
        ku=BB(bi,6)+1 ;    %upper y index of facet 1
        
        if fi==4 %north, facing north
            jt=BB(bi,4);
            y=yb(jt+1);
        else  %south, facing south
            jt=BB(bi,3);
            y=yb(jt);
        end
        
        ndim=((iu-il)+1)*((ku-kl)+1);
    end
end

%return xyz-coordinates of center and 4 corners clockwise from bottom left
%move coordinates of corners very slightly to the interior of the facet (by
%delta)
if (ci>=0 || ci<=-100) %it's not a floor -> it's a building or bounding wall
    
    if (fi==1)
        co=[(xb(iu+1)+xb(il))*.5, (yb(ju+1)+yb(jl))*.5  , z; ... %center
            xb(il)+delta        , yb(jl)+delta          , z; ... %4 corners
            xb(il)+delta        , yb(ju+1)-delta        , z; ...
            xb(iu+1)-delta      , yb(ju+1)-delta        , z; ...
            xb(iu+1)-delta      , yb(jl)+delta          , z; ...
            xb(il)        , yb(jl)          , z; ... %4 corners
            xb(il)        , yb(ju+1)        , z; ...
            xb(iu+1)      , yb(ju+1)        , z; ...
            xb(iu+1)      , yb(jl)          , z; ...
            ];
        area=(xb(iu+1)-xb(il))*(yb(ju+1)-yb(jl));
        
    elseif ((fi==2) || (fi==3))
        
        co=[ x  , (yb(ju+1)+yb(jl))*.5  , (zb(ku+1)+zb(kl))*.5; ... %center
            x  , yb(jl)+delta          , zb(kl)+delta        ; ... %4 corners
            x  , yb(jl)+delta          , zb(ku+1)-delta      ; ...
            x  , yb(ju+1)-delta        , zb(ku+1)-delta      ; ...
            x  , yb(ju+1)-delta        , zb(kl)+delta        ; ...
            x  , yb(jl)          , zb(kl)        ; ... %4 corners
            x  , yb(jl)          , zb(ku+1)      ; ...
            x  , yb(ju+1)        , zb(ku+1)      ; ...
            x  , yb(ju+1)        , zb(kl)        ; ...
            ];
        area=(zb(ku+1)-zb(kl))*(yb(ju+1)-yb(jl));
        
    elseif ((fi==4) || (fi==5))
        
        
        co=[(xb(iu+1)+xb(il))*.5     , y  , (zb(ku+1)+zb(kl))*.5  ; ... %center
            xb(il)+delta             , y  , zb(kl)+delta          ; ... %4 corners
            xb(il)+delta             , y  , zb(ku+1)-delta        ; ...
            xb(iu+1)-delta           , y  , zb(ku+1)-delta        ; ...
            xb(iu+1)-delta           , y  , zb(kl)+delta          ; ...
            xb(il)             , y  , zb(kl)          ; ... %4 corners
            xb(il)             , y  , zb(ku+1)        ; ...
            xb(iu+1)           , y  , zb(ku+1)        ; ...
            xb(iu+1)           , y  , zb(kl)          ; ...
            ];
        area=(xb(iu+1)-xb(il))*(zb(ku+1)-zb(kl));
    end
    
else  % it is a floor, not a building
    il=G(bi,1);     %lower x index of floor facet 1
    iu=G(bi,2);     %upper x index of floor facet 1
    jl=G(bi,3);     %lower y index of floor facet 1
    ju=G(bi,4);     %upper y index of floor facet 1
    z=0;
    ndim=((iu-il)+1)*((ju-jl)+1);
    co=[(xb(iu+1)+xb(il))*.5, (yb(ju+1)+yb(jl))*.5  , z; ... %center
        xb(il)+delta        , yb(jl)+delta          , z; ... %4 corners
        xb(il)+delta        , yb(ju+1)-delta        , z; ...
        xb(iu+1)-delta      , yb(ju+1)-delta        , z; ...
        xb(iu+1)-delta      , yb(jl)+delta          , z; ...
        xb(il)        , yb(jl)          , z; ... %4 corners
        xb(il)        , yb(ju+1)        , z; ...
        xb(iu+1)      , yb(ju+1)        , z; ...
        xb(iu+1)      , yb(jl)          , z; ...
        ];
    area=(xb(iu+1)-xb(il))*(yb(ju+1)-yb(jl));
    %if it is 1*1 nothing has to be done since it matches both walls
    if ju-jl>0  %floor facet along a west or east wall
        %test if it is a wall-wall-floor corner, if so, the floor facet needs
        %to be cut diagonally, see "createfloors.m"
        if cornm(iu,ju)==10  %south/west corner
            %JUST REMOVE IT ONE dx OR dy FROM THE OTHER BLOCK
            co(2+1,2)=yb(ju)+delta;
            co(6+1,2)=yb(ju);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        elseif cornm(iu,ju)==15  %south/east corner
            co(3+1,2)=yb(ju)+delta;
            co(7+1,2)=yb(ju);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        end
        if cornm(iu,jl)==8  %north/west corner
            co(1+1,2)=yb(jl+1)-delta;
            co(5+1,2)=yb(jl+1);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        elseif cornm(iu,jl)==12  %north/east corner
            co(4+1,2)=yb(jl+1)-delta;
            co(8+1,2)=yb(jl+1);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        end
        
    elseif iu-il>0 %floor facet along along a south or north wall
        if cornm(iu,ju)==10  %south/west corner
            co(4+1,1)=xb(iu)+delta;
            co(8+1,1)=xb(iu);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        elseif cornm(iu,ju)==8  %north/west corner
            co(3+1,1)=xb(iu)+delta;
            co(7+1,1)=xb(iu);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        end
        if cornm(il,ju)==12  %north/east corner
            co(2+1,1)=xb(il+1)-delta ;
            co(6+1,1)=xb(il+1) ;
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        elseif cornm(il,ju)==15  %south/east corner
            co(1+1,1)=xb(il+1)-delta;
            co(5+1,1)=xb(il+1);
            area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
        end
    end
end
