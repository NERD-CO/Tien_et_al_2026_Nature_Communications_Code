function plot_SegmentationLocalization(finame)

load(finame);

figure;
clf;
hold on;

trisurf(trim, 'FaceColor', [0 0.75 0], 'FaceAlpha', 0.5, 'LineStyle', 'none');

for tracki = 1:size(trackzeros, 2)
    if tracki == centeri
        plot3(trackzeros(1,tracki), trackzeros(2,tracki), trackzeros(3,tracki), '+k', 'MarkerSize', 10, 'LineWidth', 4);
    end
    plot3([trackzeros(1,tracki) trackzeros(1,tracki)+trackplotlen*longaxis(1)], [trackzeros(2,tracki) trackzeros(2,tracki)+trackplotlen*longaxis(2)], [trackzeros(3,tracki) trackzeros(3,tracki)+trackplotlen*longaxis(3)], '-k', 'LineWidth', 2);
end

for neui = 1:length(neupoints)
    plot3(neupoints{neui}(1), neupoints{neui}(2), neupoints{neui}(3), '.', 'Color', neucols{neui}, 'MarkerSize', 40);
end

axis square
axis equal
xlabel('+Left -Right');
ylabel('+Post -Ant');
zlabel('+Up -Down');
grid on;
box on;
xticks([5 15 25]);
yticks([0 10]);
zticks([25 35]);

set(gca,'TickLength',[0 0]);

view(245, 15);

xlabel('');
ylabel('');
zlabel('');
xticklabels([]);
yticklabels([]);
zticklabels([]);