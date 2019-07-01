
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
%pf1sf2o=single(pf1sf2);
blub=sum(vfo,2);
lblub=find(blub>1);
if ~isempty(lblub)
disp('max vf was:')
[maxvf,indexmaxvf] = max(blub)
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


ncid = netcdf.create([outputdir '/vf.nc.inp.' num2str(expnr)],'NC_WRITE');
dimidrow = netcdf.defDim(ncid,'rows',nfcts);
dimidcol = netcdf.defDim(ncid,'columns',nfcts);
varid = netcdf.defVar(ncid,'view factor','NC_FLOAT',[dimidrow dimidcol]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,vfo);
netcdf.close(ncid);

% fname = [outputdir '/vf.inp.' num2str(expnr)];
% fileID = fopen(fname,'W');
% fprintf(fileID, '# view factors between facets\n');
% fclose(fileID);
% dlmwrite(fname,vfo,'-append','delimiter',' ','precision','%4f')
% 
% ncid = netcdf.create([outputdir '/pf1sf2.nc.inp.' num2str(expnr)],'NC_WRITE');
% dimidrow = netcdf.defDim(ncid,'rows',nfcts);
% dimidcol = netcdf.defDim(ncid,'columns',nfcts);
% varid = netcdf.defVar(ncid,'percentage f1 sees of f2','NC_FLOAT',[dimidrow dimidcol]);
% netcdf.endDef(ncid);
% netcdf.putVar(ncid,varid,pf1sf2o);
% netcdf.close(ncid);

% fname = [outputdir '/pf1sf2.inp.' num2str(expnr)];
% fileID = fopen(fname,'W');
% fprintf(fileID, '# % facets see each other\n');
% fclose(fileID);
% dlmwrite(fname,pf1sf2o,'-append','delimiter',' ','precision','%4f')

fname = [outputdir '/facetarea.inp.' num2str(expnr)];
fileID = fopen(fname,'W');
fprintf(fileID, '# area of facets\n');
fclose(fileID);
dlmwrite(fname,A,'-append','delimiter',' ','precision','%4f')

end