function [mask, image ] = fillgaps( mask,image,nj,ni,pad)
%fills horizontal and vertical 1D gaps of width 1
while true
    change=false;
    for j=pad+1:nj-pad-1
        for i=pad+1:ni-pad-1
            if j==1 || i==1 %do nothing at domain edge
            continue
                
            elseif sum(sum(mask(j-1:j,i-1:i+1)))==5 && mask(j,i)==0 && mask(j+1,i)==0 %towards south. If it has 5 block neighbours, is not a block and j+1 is also not a block then fill the gap towards the south. remember that images are upside down, so j+1 is going down
                jj=j;
                while true
                    mask(jj,i)=1;
                    image(jj,i)=image(jj-1,i);
                    jj=jj+1;
                    if jj==nj %reached end of domain
                        break
                    end
                    if ~(sum(sum(mask(jj-1:jj,i-1:i+1)))==5 && mask(jj,i)==0 && mask(jj+1,i)==0)
                        break
                    end
                end
                change=true;
            elseif sum(sum(mask(j:j+1,i-1:i+1)))==5 && mask(j,i)==0 && mask(j-1,i)==0 %towards north
                jj=j;
                while true
                    mask(jj,i)=1;
                    image(jj,i)=image(jj+1,i);
                    jj=jj-1;
                    if jj==1 %reached end of domain
                        break
                    end
                    if ~(sum(sum(mask(jj:jj+1,i-1:i+1)))==5 && mask(jj,i)==0 && mask(jj-1,i)==0)
                        break
                    end
                end
                change=true;
            elseif sum(sum(mask(j-1:j+1,i-1:i)))==5 && mask(j,i)==0 && mask(j,i+1)==0 %towards east
                ii=i;
                while true
                    mask(j,ii)=1;
                    image(j,ii)=image(j,ii-1);
                    ii=ii+1;
                    if ii==ni %reached end of domain
                        break
                    end
                    if ~(sum(sum(mask(j-1:j+1,ii-1:ii)))==5 && mask(j,ii)==0 && mask(j,ii+1)==0)
                        break
                    end
                end
                change=true;
            elseif sum(sum(mask(j-1:j+1,i:i+1)))==5 && mask(j,i)==0 && mask(j,i-1)==0 %towards west
                ii=i;
                while true
                    mask(j,ii)=1;
                    image(j,ii)=image(j,ii+1);
                    ii=ii-1;
                    if ii==1 %reached end of domain
                        break
                    end
                    if ~(sum(sum(mask(j-1:j+1,ii:ii+1)))==5 && mask(j,ii)==0 && mask(j,ii-1)==0)
                        break
                    end
                end
            end
            
        end
        
        
    end
    if ~change
        break
    end
end
end


