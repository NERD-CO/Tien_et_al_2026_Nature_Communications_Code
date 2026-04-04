function plot_invim_dists(Settings)

dists2vim = Settings.Global.Locs.Segmentation.dists2vim;
postvim = Settings.Global.Locs.Segmentation.postvim;

binres = 0.2;
figure('Renderer', 'painters');
hold on;
histogram(dists2vim(dists2vim<=0), -1.6:binres:0, 'LineStyle', '-', 'LineWidth', 1, 'EdgeColor', 'k', 'FaceColor', [0 0.75 0], 'FaceAlpha', 1)
histogram(dists2vim(dists2vim>0), 0:binres:3, 'EdgeColor', 'k', 'FaceColor', 'r', 'FaceAlpha', 1, 'LineStyle', '-', 'LineWidth', 1);
histogram(dists2vim(dists2vim>0 & postvim==1), 0:binres:3, 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', 1, 'LineStyle', '-', 'LineWidth', 1);
legend({'In', 'Ant', 'Post'})

% plot([0 0], [0 16], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
nvalid = sum(~isnan(dists2vim));
xticks([-1.5 0 1.5 3]);
yticks([0 nvalid*0.05 nvalid*0.1]);
xticklabels([]);
yticklabels([]);
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.03 0.01]);