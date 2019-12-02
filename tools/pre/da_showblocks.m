function da_showblocks(bl, xh, yh, zh)

figure
view(52,23)

for i = 1:size(bl,1)     
    patch([xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1)  xh(bl(i,1))], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,4)+1) yh(bl(i,4)+1)], [zh(bl(i,6)+1)  zh(bl(i,6)+1) zh(bl(i,6)+1) zh(bl(i,6)+1)], [245 245 245] ./ 255) % [245 245 245] ./ 255
    patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,3)) yh(bl(i,3))], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
    patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,4)+1)  yh(bl(i,4)+1) yh(bl(i,4)+1) yh(bl(i,4)+1)], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
    patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,1)) ], [yh(bl(i,4)+1)  yh(bl(i,4)+1) yh(bl(i,3)) yh(bl(i,3))], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
    patch([xh(bl(i,2)+1) xh(bl(i,2)+1) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,4)+1) yh(bl(i,4)+1)], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
end

zlim([0 zh(end)]); %/(r.blockheight-1))
xlim([0 xh(end)]); %/(r.blockheight-1))
ylim([0 yh(end)]); %/(r.blockheight-1))

set(gca,'ticklabelinterpreter','latex')
xlabel('x [m]','interpreter','latex')
ylabel('y [m]','interpreter','latex')
zlabel('z [m]','interpreter','latex')
set(gca,'BoxStyle','full','Box','on')
daspect([1 1 1])

end