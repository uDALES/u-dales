function da_showfacets(fctl, xb, yb, zb, il, iu, jl, ju, kl, ku)

cmap = colormap('parula');
figure
    for i=1:nfcts
        switch fctl(i, 1)
            case {top, bot}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl))];
            case {west, east}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+kl))];
            case {north, south}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku))];
        end
        
        ci = min(floor(double(fctl(i, 2))/double(max(fctl(:, 2)))*length(cmap))+1, length(cmap));
        patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        d = [0, 0, 0]; a = 0.25;
        switch(fctl(i, 1))
            case top
                d(3) = a * dz(1);
            case west
                d(1) = -a*dx(1);
            case east
                d(1) = a*dx(1);
            case south
                d(2) = -a*dy(1);
            case north
                d(2) = a*dy(1);
        end
        text(mean(x)+d(1), mean(y)+d(2), mean(z)+d(3), num2str(i), ...
            'horizontalalignment', 'center')
        hold on
        title('Facets')
    end
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    xlim([0 xb(end)])
    ylim([0 yb(end)])
    zlim([0 zb(end)])
    
end