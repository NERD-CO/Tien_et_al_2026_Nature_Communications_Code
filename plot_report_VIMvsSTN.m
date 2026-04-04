
function varargout = plot_report_VIMvsSTN(Settings, Realtime, RealtimeSig, STN_Settings, STN_Realtime, STN_RealtimeSig);

% Doing peak speed only
aistr = 'peak';

%% Peak speed aligned modulation, VIM
npt = size(Realtime(1).(aistr).rates,1);
nall = Settings.Global.allnneu;
nalpha = length(Settings.alph2do);

allspeeds = [];
allsig = nan(nall,nalpha);
allsigofem = nan(nall, npt, nalpha);

ticker = 1;
for dati = 1:length(Realtime)

    nneu = size(Realtime(dati).(aistr).rates,3);

    allspeeds = [allspeeds; Realtime(dati).(aistr).speeds'];

    for alphi = 1:nalpha
        alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

        allsig(ticker:ticker+nneu-1,alphi) = RealtimeSig(dati).(aistr).(alphstr).sigi;
        allsigofem(ticker:ticker+nneu-1,:,alphi) = RealtimeSig(dati).(aistr).(alphstr).sigofem;
    end

    ticker = ticker+nneu;
end

meanspeed = nanmean(allspeeds,1);
stdspeed = nanstd(allspeeds,[],1);

plotx = 1:npt;
errdown = (meanspeed - stdspeed);
errup = (meanspeed + stdspeed);

%% Peak speed aligned modulation, STN
STN_npt = size(STN_Realtime(1).(aistr).rates,1);
STN_nall = STN_Settings.Global.allnneu;
STN_nalpha = length(STN_Settings.alph2do);

STN_allspeeds = [];
STN_allsig = nan(STN_nall,STN_nalpha);
STN_allsigofem = nan(STN_nall, STN_npt, STN_nalpha);

ticker = 1;
for dati = 1:length(STN_Realtime)

    nneu = size(STN_Realtime(dati).(aistr).rates,3);

    STN_allspeeds = [STN_allspeeds; STN_Realtime(dati).(aistr).speeds'];

    for alphi = 1:nalpha
        alphstr = ['alpha' num2str(100*STN_Settings.alph2do(alphi),'%0.2i')];

        STN_allsig(ticker:ticker+nneu-1,alphi) = STN_RealtimeSig(dati).(aistr).(alphstr).sigi;
        STN_allsigofem(ticker:ticker+nneu-1,:,alphi) = STN_RealtimeSig(dati).(aistr).(alphstr).sigofem;
    end

    ticker = ticker+nneu;
end

STN_meanspeed = nanmean(STN_allspeeds,1);
STN_stdspeed = nanstd(STN_allspeeds,[],1);

STN_plotx = 1:STN_npt;
STN_errdown = (STN_meanspeed - STN_stdspeed);
STN_errup = (STN_meanspeed + STN_stdspeed);



alphi = 2;
alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
%% Mini Stacked plots

figure('Renderer', 'painters'); hold on;
posit = get(gcf, 'Position');
posit(3) = 500;
set(gcf, 'Position', posit);

% Speed
sp1 = subplot(2,1,1); hold on;
title(aistr)
patch([plotx, fliplr(plotx)], [errdown, fliplr(errup)], 'k', 'FaceAlpha', 0.25, 'LineStyle', 'none'); % VIM
plot(meanspeed, '-k', 'LineWidth', 4); % VIM
patch([STN_plotx, fliplr(STN_plotx)], [STN_errdown, fliplr(STN_errup)], 'm', 'FaceAlpha', 0.25, 'LineStyle', 'none'); % STN
plot(STN_meanspeed, '-m', 'LineWidth', 4); % STN
xticks([]);
ylim([0 520]);
set(gca, 'TickDir', 'out');
ax = gca;
ax.FontSize = 30;
xlim([0.5, Settings.Realtime.naround+0.5]);
yticks([0 250 500]);
yticklabels([]);

plot([5 55], [150 150], '-k', 'LineWidth', 8); % time scale
set(gca,'TickLength',[0.03 0.01]);
% Sigs
sp2 = subplot(2,1,2); hold on;
plot(sum(allsigofem(:,:,alphi),1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'black'); % VIM
plot(sum(STN_allsigofem(:,:,alphi),1)./STN_nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'm'); % STN

xticks([0.5 0.5+Settings.Realtime.naround/2 0.5+Settings.Realtime.naround])
set(gca, 'XTickLabels', []);

a = get(gca, 'XTickLabels');
set(gca, 'XTickLabels', a, 'FontSize', 30);
set(gca, 'TickDir', 'out');

xlim([0.5, Settings.Realtime.naround+0.5]);

if alphi == 1
    ylim([0 0.35]);
    yticks([0 .1 .2 .3]);
else
    ylim([0 0.25]);
    yticks([0 0.1 0.2]);
end
yticklabels([]);
set(gca,'TickLength',[0.03 0.01]);

%% Now do the percentiles and plot
alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
timevec = -0.495:0.01:0.495;

%% VIM
vimsum = squeeze(sum(allsigofem(:,:,alphi),1));
vimcdf = cumsum(vimsum./nall);
vimtimes = [];
for vi = 1:length(vimsum)
    vimtimes = [vimtimes; ones(vimsum(vi),1)*timevec(vi)];
end
mvim = median(vimtimes);

vimanysig = sum(any(allsigofem(:,:,alphi),2));
vimnall = nall;
vimpctsig = vimanysig/vimnall;

%% STN
stnsum = squeeze(sum(STN_allsigofem(:,:,alphi),1));
stncdf = cumsum(stnsum./nall);
stntimes = [];
for si = 1:length(stnsum)
    stntimes = [stntimes; ones(stnsum(si),1)*timevec(si)];
end
mstn = median(stntimes);

stnanysig = sum(any(STN_allsigofem(:,:,alphi),2));
stnnall = nall;
stnpctsig = stnanysig/stnnall;

figure('Renderer', 'painters'); hold on;
posit = get(gcf, 'Position');
posit(3) = 500;
posit = posit*1.5;
set(gcf, 'Position', posit);
axis square
hold on;

% Plot diagonal
plot([timevec(1) timevec(end)], [0 1], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 2);
plot([0 0], [0 1], '--', 'Color', [0.75 0.75 0.75]);

% Plot medians
% plot([0 mvim], [0.5 0.5], '--k');
plot([mvim mvim], [0 0.5], '--k');
% plot([0 mstn], [0.5 0.5], '--m');
plot([mstn mstn], [0 0.5], '--m');

% Plot stairs
plot([timevec(1) timevec(end)], [0.5 0.5], '--', 'Color', [0.75 0.75 0.75]);
stairs(timevec(1):Settings.Realtime.step:timevec(end), vimcdf/max(vimcdf), '-k', 'LineWidth', 2);
plot([0 0], [0 vimcdf(1)/max(vimcdf)], '-k', 'LineWidth', 2);
stairs(timevec(1):Settings.Realtime.step:timevec(end), stncdf/max(stncdf), '-m', 'LineWidth', 2);
plot([0 0], [0 stncdf(1)/max(stncdf)], '-m', 'LineWidth', 2);


lims = [-Settings.Realtime.naround*Settings.Realtime.step/2, Settings.Realtime.naround*Settings.Realtime.step/2];
xticks([lims(1) mstn 0 mvim lims(end)]);
xlim([lims(1) lims(end)]);
yticks([0 0.5 1]);
ylim([0 1]);
xticklabels([]);
yticklabels([]);
set(gca,'TickLength',[0.03 0.01]);
set(gca, 'TickDir', 'out');

% Return values to report
R.vimanysig = vimanysig;
R.vimnall = vimnall;
R.vimpctsig = vimpctsig;
R.stnanysig = stnanysig;
R.stnnall = stnnall;
R.stnpctsig = stnpctsig;
R.mvim = mvim;
R.mstn = mstn;
R.pVSRankSum = ranksum(stntimes, vimtimes);
R.pVIMSignRank = signrank(vimtimes);
R.pSTNSignRank = signrank(stntimes);

varargout{1} = R;