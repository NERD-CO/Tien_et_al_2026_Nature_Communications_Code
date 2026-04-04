function plot_report_CompareInVim(finame)

load(finame)

fn = fieldnames(comparevim);
for fni = 1:length(fn)
    eval([fn{fni} ' = comparevim.' fn{fni} ';']);
end

figure; clf; hold on;
posit = get(gcf, 'Position');
% set(gcf, 'Position', posit/2);

plot([timevec(1) timevec(end)], [0 1], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 2);
plot([0.4 0.4], [0 1], '--', 'Color', [0.75 0.75 0.75]);
plot([timevec(1) timevec(end)], [0.5 0.5], '--', 'Color', [0.75 0.75 0.75]);
plot([inmed inmed], [0 0.5], '--', 'Color', [0 0.75 0]);
plot([antmed antmed], [0 0.5], '--', 'Color', 'r');
plot([postmed postmed], [0 0.5], '--', 'Color', 'b');
stairs(timevec, cdfs(:,1), '-', 'Color', [0 0.75 0], 'LineWidth', 2);
stairs(timevec, cdfs(:,2), '-r', 'LineWidth', 2);
stairs(timevec, cdfs(:,3), '-b', 'LineWidth', 2); 
xticks([-.4 0 0.8 1.2]);
xticklabels([]);
yticks([0 0.5 1]);
yticklabels([]);
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.03 0.01]);
axis square;

figure;
incols = {[0 0.75 0], 'r', 'b'};

subplot(1,3,1); hold on;
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.03 0.01]);
% totprop = sum(nsig_invim)/sum(nall_invim);
alphi = 2;
for typei = 1:nclass
    bar(typei, sigprops(typei,alphi), 'EdgeColor', 'k', 'FaceColor', incols{typei})
    text(typei, 0, ['n=' num2str(nall_invim(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'w');
end
% plot([0 nclass+1], [totprop totprop], '--r');
xlim([0 nclass+1]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xticks(1:nclass);
xticklabels({'InVim', 'Ant', 'Post'});
% text(nclass+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'r');
if isempty(pairsigschi{alphi})
    text(mean([1 nclass]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
else
    for pairi = 1:size(pairsigschi{alphi},1)
        plot([pairsigschi{alphi}(1) pairsigschi{alphi}(2)], [1 1]-pairi*0.05, '-k');
    end
end

subplot(1,3,2); hold on;
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.03 0.01]);
% dirtotprop = sum(ndirsig_invim)/sum(nall_invim);
alphi = 1;
for typei = 1:nclass
    bar(typei, dirsigprops(typei,alphi), 'EdgeColor', 'k', 'FaceColor', incols{typei})
    text(typei, 0, ['n=' num2str(nall_invim(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'w');
end
% plot([0 nclass+1], [totprop totprop], '--r');
xlim([0 nclass+1]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xticks(1:nclass);
xticklabels({'InVim', 'Ant', 'Post'});
% text(nclass+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'r');
if isempty(dirpairsigschi{alphi})
    text(mean([1 nclass]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
else
    for pairi = 1:size(dirpairsigschi{alphi},1)
        plot([dirpairsigschi{alphi}(1) dirpairsigschi{alphi}(2)], [1 1]-pairi*0.05, '-k');
    end
end

subplot(1,3,3); hold on;
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.03 0.01]);
% lagtotprop = sum(nlagsig_invim)/sum(nall_invim);
alphi = 2;
for typei = 1:nclass
    bar(typei, lagsigprops(typei, alphi), 'EdgeColor', 'k', 'FaceColor', incols{typei})
    text(typei, 0, ['n=' num2str(nall_invim(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'w');
end
% plot([0 nclass+1], [totprop totprop], '--r');
xlim([0 nclass+1]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xticks(1:nclass);
xticklabels({'InVim', 'Ant', 'Post'});
% text(nclass+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'r');
if isempty(lagpairsigschi{alphi})
    text(mean([1 nclass]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
else
    for pairi = 1:size(lagpairsigschi{alphi},1)
        plot([lagpairsigschi{alphi}(1) lagpairsigschi{alphi}(2)], [1 1]-pairi*0.05, '-k');
    end
end
