% clear all
% close all

%% can facet 1 see facet 2

centerweight=0.5;  %if center is sunlit, centerweight% of the facet are sunlit
cornerweight=(1-centerweight)/4; %if 1 corner is sunlit, cornerweight% of the facet is sunlit

%% read files

% CHECK NHEADER NUMBERS WITH TEXT FILES
%blocks
nheader=2;
try %in case file is empty -> no blocks
BB = dlmread([outputdir '/blocks.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
catch
BB =[];
end

%blocks for intersection
nheader=3;
try %in case file is empty -> no blocks
B = dlmread([tempdir '/bbri.inp'],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
catch
B =[];
end
%floors
nheader=3;
%G = dlmread(['floors.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju  type
G = dlmread([outputdir '/floors.txt'],'',nheader,0);  %#il   iu   jl    ju  type
%bounding walls
nheader=3;
W = dlmread([outputdir '/boundingwalls.txt'],'',nheader,0);  %#il   iu   jl    ju  type
%facets
nheader=1;
F = dlmread([outputdir '/facets.inp.' num2str(expnr)],'',nheader,0);  %#   or     wl    blk    bld


[nblocks, nbi]=size(B);
[nfcts, nfacprop]=size(F);
[nfl, nflprop]=size(G);
[nbw, nbwprop]=size(W);

%% create blocks to test for intersection
%
bl=zeros(nblocks+nbw,6);
for k=1:nblocks
    xl=xb(B(k,1));
    xu=xb(B(k,2)+1);
    yl=yb(B(k,3));
    yu=yb(B(k,4)+1);
    zl=zb(B(k,5)+1);
    zu=zb(B(k,6)+2);
    bl(k,:)=[xl xu yl yu zl zu];
end
for k=1:nbw
    fi3=F(k+(nfcts-nfl-nbw),1);
    il=W(k,1);
    iu=W(k,2);
    jl=W(k,3);
    ju=W(k,4);
    kl=W(k,5)+1;
    ku=W(k,6)+1;
    if (fi3==2)
        xl=xb(end);
        xu=xb(end)+0.1;
        yl=yb(jl);
        yu=yb(ju+1);
        zl=zb(kl);
        zu=zb(ku+1);
    elseif (fi3==3)
        xl=xb(1)-0.1;
        xu=xb(1);
        yl=yb(jl);
        yu=yb(ju+1);
        zl=zb(kl);
        zu=zb(ku+1);
    elseif (fi3==4)
        xl=xb(il);
        xu=xb(iu+1);
        yl=yb(1)-0.1;
        yu=yb(1);
        zl=zb(kl);
        zu=zb(ku+1);
    else %if (fi==5)
        xl=xb(il);
        xu=xb(iu+1);
        yl=yb(end);
        yu=yb(end)+0.1;
        zl=zb(kl);
        zu=zb(ku+1);
    end
    bl(k+nblocks,:)=[xl xu yl yu zl zu];
end
gl=zeros(nfl,6); %[xl xu yl yu zl zu] %space coordinates of floors
for j=1:nfl
    xl=xb(G(j,1));
    xu=xb(G(j,2)+1);
    yl=yb(G(j,3));
    yu=yb(G(j,4)+1);
    zl=0-0.1;
    zu=0;
    gl(j,:)=[xl xu yl yu zl zu];
end





%% View factors
% Calculate fraction of facets seeing each other (not blocked)
% Calculate orientation, no overlap on orthogonal planes
% calculate viewfactors

vf=zeros(nfcts,nfcts);

pf1sf2=zeros(nfcts,nfcts);
A=zeros(nfcts,1);
tim=zeros(nfcts-1,1);
%%
tic
%do dummy parallelisation by starting multiple matlab and only iterate over
%a fraction and then write to file and merge files
for i=1:(nfcts-1)
    bi=F(i,3); %block index
    fi=F(i,1); %facet index
    ci=F(i,4); %building index (-1 for roads, -99 for bounding wall)
    [ ndima, areaa, coa] = detsub(i,F,BB,G,W,cornm,xb,yb,zb,delta);
    for j=(i+1):nfcts
        bi2=F(j,3); %block index
        fi2=F(j,1); %facet index
        ci2=F(j,4); %building index (-1 for roads, -99 for bounding wall)
        if ((fi2==fi) || ((ci2*bi2)==(ci*bi)))  %facets looking in same direction or on same block&building (ci*bi is unique), bi alone could also be a floor or a bounding wall
            continue
        end
        
        [ ndimb, areab, cob] = detsub(j,F,BB,G,W,cornm,xb,yb,zb,delta);
        [ prblckdij ] = prblckd(i,j,coa,ndima,false,-999,cob,ndimb,F,centerweight,cornerweight,nblocks,nbw,bl);
        c=1-prblckdij;
        pf1sf2(i,j)=c;
        pf1sf2(j,i)=c;
        A(i)=areaa;
        A(j)=areab;
        
        if prblckdij<0.98 %view is not blocked, vf calclation
            if (prblckdij>=(2*cornerweight-eps)) %at least two corners have to be blocked, otherwise it cannot possibly be blocked by itself
                %slice up
                %pf1sf2u also has to include check for center, in case center is
                %blocked
                [coaa, cobb, pf1sf2u] = slice(fi,fi2,coa,cob,cornerweight,centerweight);
                c=c+pf1sf2u;
                pf1sf2(i,j)=c;
                pf1sf2(j,i)=c;
            else
                coaa=coa(6:9,:);
                cobb=cob(6:9,:);
            end
            
            %determine if they form a corner/have common edge
            %determine the vertices pair that forms that corner
            vcorner=[0 0 1]; %corner vertices for i and j & sign of correction (as always, clockwise from bottom left)
            glpo=6; %order of gauss-legendre polynomial
            glpop=15; %more accurate, order of gauss-legendre polynomial
            if (fi==1) %fi is horizontal, test if corner with a vertical wall
                if  (fi2==2)
                    if all(coaa(4,:)==cobb(1,:)) && all(coaa(3,:)==cobb(4,:)) %they share one edge completely
                        vcorner=[3 4 -1];
                        glpo=glpop;
                        % elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,3)==cobb(1,3)) %they potentially share parts of an edge
                        %     glpo=18;  %NEEDS TO BE AN EVEN NUMBER!!!
                    end
                elseif (fi2==3)
                    if all(coaa(1,:)==cobb(1,:)) && all(coaa(2,:)==cobb(4,:)) %they share one edge completely
                        vcorner=[1 4 1];
                        glpo=glpop;
                        % elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                        %     glpo=18;
                    end
                elseif (fi2==4)
                    if all(coaa(1,:)==cobb(1,:)) && all(coaa(4,:)==cobb(4,:)) %they share one edge completely
                        corner=[4 4 -1];
                        glpo=glpop;
                        %  elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                        %      glpo=18;
                    end
                elseif (fi2==5)
                    if all(coaa(2,:)==cobb(1,:)) && all(coaa(3,:)==cobb(4,:)) %they share one edge completely
                        vcorner=[2 4 1];
                        glpo=glpop;
                        %  elseif (coaa(2,2)==cobb(1,2)) && (coaa(2,3)==cobb(1,3)) %they potentially share parts of an edge
                        %      glpo=18;
                    end
                end
            elseif (fi==2)
                if  (fi2==4)
                    if all(coaa(1,:)==cobb(4,:)) && all(coaa(2,:)==cobb(3,:)) %they share one edge completely
                        vcorner=[1 3 1];
                        glpo=glpop;
                     elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,2)==cobb(4,2)) %they potentially share parts of an edge
                         1
                         glpo=18;
                    end
                elseif (fi2==5)
                    if all(coaa(4,:)==cobb(4,:)) && all(coaa(3,:)==cobb(3,:)) %they share one edge completely
                        vcorner=[3 3 -1];
                        glpo=glpop;
                     elseif (coaa(4,1)==cobb(4,1)) && (coaa(4,2)==cobb(4,2)) %they potentially share parts of an edge
                        2
                         glpo=18;
                    end
                elseif (fi2==1)
                    if all(coaa(1,:)==cobb(4,:)) && all(coaa(4,:)==cobb(3,:)) %they share one edge completely
                        vcorner=[4 3 -1];
                        glpo=glpop;
                        %                    elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,3)==cobb(4,3)) %they potentially share parts of an edge
                        %                        glpo=18;
                    end
                end
            elseif (fi==3)
                if  (fi2==4)
                    if all(coaa(1,:)==cobb(1,:)) && all(coaa(2,:)==cobb(2,:)) %they share one edge completely
                        vcorner=[1 1 -1];
                        glpo=glpop;
                     elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,2)==cobb(1,2)) %they potentially share parts of an edge
                         3
                         glpo=18;
                    end
                elseif (fi2==5)
                    if all(coaa(4,:)==cobb(1,:)) && all(coaa(3,:)==cobb(2,:)) %they share one edge completely
                        vcorner=[3 1 1];
                        glpo=glpop;
                     elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,2)==cobb(1,2)) %they potentially share parts of an edge
                         4
                         glpo=18;
                    end
                elseif (fi2==1)
                    if all(coaa(1,:)==cobb(1,:)) && all(coaa(4,:)==cobb(2,:)) %they share one edge completely
                        vcorner=[4 1 1];
                        glpo=glpop;
                        %  elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                        %      glpo=18;
                    end
                end
            elseif (fi==4)
                if  (fi2==2)
                    if all(coaa(4,:)==cobb(1,:)) && all(coaa(3,:)==cobb(2,:)) %they share one edge completely
                        vcorner=[3 1 1];
                        glpo=glpop;
                     elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,2)==cobb(1,2)) %they potentially share parts of an edge
                         5
                         glpo=18;
                    end
                elseif (fi2==3)
                    if all(coaa(1,:)==cobb(1,:)) && all(coaa(2,:)==cobb(2,:)) %they share one edge completely
                        vcorner=[1 1 -1];
                        glpo=glpop;
                     elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,2)==cobb(1,2)) %they potentially share parts of an edge
                         6
                         glpo=18;
                    end
                elseif (fi2==1)
                    if all(coaa(1,:)==cobb(1,:)) && all(coaa(4,:)==cobb(4,:)) %they share one edge completely
                        vcorner=[4 4 -1];
                        glpo=glpop;
                        %    elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                        %       glpo=18;
                    end
                end
            elseif (fi==5)
                if  (fi2==2)
                    if all(coaa(4,:)==cobb(4,:)) && all(coaa(3,:)==cobb(3,:)) %they share one edge completely
                        vcorner=[3 3 -1];
                        glpo=glpop;
                     elseif (coaa(4,1)==cobb(4,1)) && (coaa(4,2)==cobb(4,2)) %they potentially share parts of an edge
                         7
                         glpo=18;
                    end
                elseif (fi2==3)
                    if all(coaa(1,:)==cobb(4,:)) && all(coaa(2,:)==cobb(3,:)) %they share one edge completely
                        vcorner=[1 3 1];
                        glpo=glpop;
                     elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,2)==cobb(4,2)) %they potentially share parts of an edge
                         8
                         glpo=18;
                    end
                elseif (fi2==1)
                    if all(coaa(1,:)==cobb(2,:)) && all(coaa(4,:)==cobb(3,:)) %they share one edge completely
                        vcorner=[4 2 1];
                        glpo=glpop;
                        % elseif (coaa(1,2)==cobb(2,2)) && (coaa(1,3)==cobb(2,3)) %they potentially share parts of an edge
                        %     glpo=18
                    end
                end
            end
            %             if (fi==1) %fi is horizontal, test if corner with a vertical wall
            %              %e.g   j is west wall,      ixu==jx                    iz==jzl                    iyl>=jyl              iyl<=jyu                  iyu>=jyl              iyu<=jyu                  iyu>=jyl              iyl<=jyl                  iyu>=jyu              iyl<=jyu
            %                 if  ( (fi2==2) && (coaa(3,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) && ( (coaa(1,2)>=cobb(1,2)&&coaa(1,2)<=cobb(4,2)) || (coaa(2,2)>=cobb(1,2)&&coaa(2,2)<=cobb(4,2)) || (coaa(2,2)>=cobb(1,2)&&coaa(1,2)<=cobb(1,2)) || (coaa(2,2)>=cobb(4,2)&&coaa(1,2)<=cobb(4,2))) )
            %                     vcorner=[3 4 -1];
            %                 elseif ((fi2==3) && (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) && ( (coaa(1,2)>=cobb(1,2)&&coaa(1,2)<=cobb(4,2)) || (coaa(2,2)>=cobb(1,2)&&coaa(2,2)<=cobb(4,2)) || (coaa(2,2)>=cobb(1,2)&&coaa(1,2)<=cobb(1,2)) || (coaa(2,2)>=cobb(4,2)&&coaa(1,2)<=cobb(4,2))) )
            %                     vcorner=[1 4 1];
            %                 elseif ((fi2==4) && (coaa(1,2)==cobb(1,2)))
            %                     vcorner=[4 4 -1];
            %                 elseif ((fi2==5) && (coaa(2,2)==cobb(1,2)))
            %                     vcorner=[2 4 1];
            %                 end
            %             elseif (fi==2)
            %                 if     ((fi2==4) && (coaa(1,2)==cobb(1,2)))
            %                     vcorner=[1 3 1];
            %                 elseif ((fi2==5) && (coaa(3,2)==cobb(1,2)))
            %                     vcorner=[3 3 -1];
            %                 elseif ((fi2==1) && (coaa(1,3)==cobb(1,3)))
            %                     vcorner=[4 3 1];
            %                 end
            %             elseif (fi==3)
            %                 if     ((fi2==4) && (coaa(1,2)==cobb(1,2)))
            %                     vcorner=[1 1 -1];
            %                 elseif ((fi2==5) && (coaa(3,2)==cobb(1,2)))
            %                     vcorner=[3 1 1];
            %                 elseif ((fi2==1) && (coaa(1,3)==cobb(1,3)))
            %                     vcorner=[4 1 1];
            %                 end
            %             elseif (fi==4)
            %                 if     ((fi2==2) && (coaa(3,1)==cobb(1,1)))
            %                     vcorner=[3 1 1];
            %                 elseif ((fi2==3) && (coaa(1,1)==cobb(1,1)))
            %                     vcorner=[1 1 -1];
            %                 elseif ((fi2==1) && (coaa(1,3)==cobb(1,3)))
            %                     vcorner=[4 4 -1];
            %                 end
            %             elseif (fi==5)
            %                 if     ((fi2==2) && (coaa(3,1)==cobb(1,1)))
            %                     vcorner=[3 3 -1];
            %                 elseif ((fi2==3) && (coaa(1,1)==cobb(1,1)))
            %                     vcorner=[1 3 1];
            %                 elseif ((fi2==1) && (coaa(1,3)==cobb(1,3)))
            %                     vcorner=[4 2 1];
            %                 end
            %             end
            
            
            %calculate viewfactor
            [F12,F21]=ViewFactor(coaa,cobb,areaa,areab,glpo,vcorner);
            
            %increase glpo if view factors are big
            if F12>0.5 || F21>0.5
                [F12,F21]=ViewFactor(coaa,cobb,areaa,areab,glpo+10,vcorner);
            end
            
            
            % round to 1promille accuracy? i.e. cut smaller 1promille completely
            % this violates reciprocity though
            
             vf(i,j)=F12*c;
             vf(j,i)=F21*c;
         %   vf(i,j)=floor(F12*c*100)/100;
          %  vf(j,i)=floor(F21*c*100)/100;
        end
    end
    i
    toc
    tim(i)=toc;
end
disp('done calculating vf')
%%
vfo=single(vf);
pf1sf2o=single(pf1sf2);
blub=sum(vfo,2);
lblub=find(blub>1);
if ~isempty(lblub)
disp('max vf was:')
max(blub)
for i=1:length(lblub)
vfo(lblub(i),:)=vfo(lblub(i),:)/blub(lblub(i));
end
end
svf=max(1-sum(vfo,2),0);
% end
%% write
if lwritefiles
% save([outputdir '/vfo'],'vfo')
% save([outputdir '/pf1sf2o'],'pf1sf2o')
% save([outputdir '/A'],'A')

fname = [outputdir '/svf.inp.' num2str(expnr)];
fileID = fopen(fname,'W');
fprintf(fileID, '# sky view factors\n');
fclose(fileID);
dlmwrite(fname,svf,'-append','delimiter',' ','precision','%4f')


fname = [outputdir '/vf.inp.' num2str(expnr)];
fileID = fopen(fname,'W');
fprintf(fileID, '# view factors between facets\n');
fclose(fileID);
dlmwrite(fname,vfo,'-append','delimiter',' ','precision','%4f')

fname = [outputdir '/pf1sf2.inp.' num2str(expnr)];
fileID = fopen(fname,'W');
fprintf(fileID, '# % facets see each other\n');
fclose(fileID);
dlmwrite(fname,pf1sf2o,'-append','delimiter',' ','precision','%4f')

fname = [outputdir '/facetarea.inp.' num2str(expnr)];
fileID = fopen(fname,'W');
fprintf(fileID, '# area of facets\n');
fclose(fileID);
dlmwrite(fname,A,'-append','delimiter',' ','precision','%4f')
%
%
% fname = ['tvf.inp.' num2str(expnr)];
% fileID = fopen(fname,'w');
% fprintf(fileID, 'environmental and sky view factors\n');
% fclose(fileID);
% dlmwrite(fname,tvf,'-append','delimiter',' ','precision','%4f')
%
% fname = 'facet-area.txt';
% fileID = fopen(fname,'w');
% fprintf(fileID, 'facet areas [m2]\n');
% fclose(fileID);
% dlmwrite(fname ,area,'-append','delimiter',' ','precision','%4f')
end