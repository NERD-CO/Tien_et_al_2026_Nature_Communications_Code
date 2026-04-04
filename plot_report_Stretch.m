function varargout = plot_report_Stretch(Settings, Stretch, StretchSig)

nsess = length(Stretch);

npt = Settings.Stretch.nbef + Settings.Global.nbet + Settings.Stretch.naft;
nall = Settings.Global.allnneu;
nalpha = length(Settings.alph2do);

allmrates = nan(nall, npt);
allnormrates = nan(nall, npt);
allspeeds = [];
allsig = nan(nall,nalpha);
allallofem = nan(nall, npt, nalpha);
allsigofem = nan(nall, npt, nalpha);
allsighilo = nan(nall,npt,nalpha);

ticker = 1;
for si = 1:nsess
    nneu = size(Stretch(si).rates,3);

    allmrates(ticker:ticker+nneu-1,:) = squeeze(mean(Stretch(si).rates,2))';
    allspeeds = [allspeeds; Stretch(si).speeds'];

    for alphi = 1:nalpha
        alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

        allsig(ticker:ticker+nneu-1,alphi) = StretchSig(si).(alphstr).sigi;
        allallofem(ticker:ticker+nneu-1,:,alphi) = StretchSig(si).(alphstr).ps < Settings.alph2do(alphi);
        allsigofem(ticker:ticker+nneu-1,:,alphi) = StretchSig(si).(alphstr).sigofem;
        allsighilo(ticker:ticker+nneu-1,:,alphi) = StretchSig(si).(alphstr).sighilo;
    end

    ticker = ticker+nneu;
end

allspeeds = allspeeds(~any(allspeeds>800,2),:);
meanspeed = nanmean(allspeeds,1);
stdspeed = nanstd(allspeeds,[],1);

plotx = 1:npt;
errdown = (meanspeed - stdspeed);
errup = (meanspeed + stdspeed);

for alli = 1:nall
    allnormrates(alli,:) = (allmrates(alli,:) - min(allmrates(alli,:)))/(range(allmrates(alli,:)));
end

% Peri-reach windows
periwin = [false(1,Settings.Stretch.nbef-Settings.Stretch.periwin) true(1,Settings.Stretch.periwin) true(1,Settings.Global.nbet) true(1,Settings.Stretch.periwin) false(1,Settings.Stretch.naft-Settings.Stretch.periwin)];
% These windows are for AFTER perireach window has been applied. Pre 400ms, first
% half of reach, 2nd half of reach, post 400ms
segwin = [true(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) true(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) true(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) true(1,Settings.Stretch.periwin)];

alphi = 2;
alphstr = ['alpha' num2str(100*Settings.alph2do(alphi), '%0.2i')];
sigi = any(allsigofem(:,periwin,alphi),2);
nsig = sum(sigi);
normrates = allnormrates(sigi,:);
sigofem = allsigofem(sigi,:,alphi);
sighilo = allsighilo(sigi,:,alphi);
[~,peakt] = max(normrates,[],2);
[~,sorti] = sort(peakt);
sorti = flipud(sorti);
normrates = normrates(sorti,:);
sigofem = sigofem(sorti,:);
sighilo = sighilo(sorti,:);

% Do some stats!
winsigofem = sigofem(:,periwin);
winsighilo = sighilo(:,periwin);

% Count sig samples for chi square test
totcount = sum(sum(winsigofem));
segcount = nan(4,1);
segexpected = nan(4,1);
segunits = nan(4,1);
for segi = 1:4
    segcount(segi) = sum(sum(winsigofem(:,segwin(segi,:))));
    segexpected(segi) = totcount*(sum(segwin(segi,:))/size(segwin,2));
    segunits(segi) = sum(any(winsigofem(:,segwin(segi,:)),2));
end

inunits = sum(any(winsigofem(:,segwin(2,:)),2) | any(winsigofem(:,segwin(3,:)),2));

% Do chisq test for in reach vs out of reach
incount = sum(segcount(2:3));
outcount = sum(segcount([1 4]));
inexpected = sum(segexpected(2:3));
outexpected = sum(segexpected([1 4]));
inoutchi = (((incount-inexpected).^2)./inexpected) + (((outcount-outexpected).^2)./outexpected);
inoutchip = chi2cdf(inoutchi,1,'upper');

% Pairwise chisq tests
npairs = nchoosek(4,2);
pairchip = nan(npairs,1);
segpairs = nan(npairs,2);
ticker = 1;
for segi1 = 1:3
    for segi2 = (segi1+1):4
        count1 = segcount(segi1);
        count2 = segcount(segi2);
        expected1 = (count1+count2)*(sum(segwin(segi1,:)))/(sum(segwin(segi1,:))+sum(segwin(segi2,:)));
        expected2 = (count1+count2)*(sum(segwin(segi2,:)))/(sum(segwin(segi1,:))+sum(segwin(segi2,:)));
        pairchi = (((count1-expected1).^2)./expected1) + (((count2-expected2).^2)./expected2);
        pairchip(ticker) = chi2cdf(pairchi,1,'upper');
        segpairs(ticker,:) = [segi1 segi2];
        ticker = ticker+1;
    end
end

R(alphi).nsig = nsig;
R(alphi).sigi = sigi;
R(alphi).nsigsamples = totcount;
R(alphi).npossamples = sum(sum(winsighilo==1));
R(alphi).onlypos = all(winsighilo >= 0, 2);
R(alphi).onlyneg = all(winsighilo <= 0, 2);
R(alphi).posneg = any(winsighilo == -1, 2) & any(winsighilo == 1, 2);
R(alphi).segcount = segcount;
R(alphi).incount = incount;
R(alphi).inoutchip = inoutchip;
R(alphi).pairchip = pairchip;
R(alphi).segpairs = segpairs;
R(alphi).winsigofem = winsigofem;
R(alphi).segwin = segwin;
R(alphi).segunits = segunits;
R(alphi).inunits = inunits;
[R(alphi).npeakmod, R(alphi).peakmodtime] = max(squeeze(sum(sigofem,1)));

%% Make the fancy subplots

% FR Bars plot
figure('Position', [10 10 700 1000]);
clf;
hold on;

colormap gray

imalpha = ones(size(normrates));
imalpha(isnan(normrates)) = 0;
imagesc(normrates, 'AlphaData', imalpha);
xticks([0.5, Settings.Stretch.nbef+0.5-Settings.Stretch.periwin, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5+Settings.Stretch.periwin, npt+0.5])
xticklabels({});

xlim([0.5 plotx(end)+0.5]);
ylim([0.5 nsig+0.5]);
ax = gca;
ax.FontSize = 30;
set(gca, 'TickDir', 'out');
yticks([]);
yticklabels([]);
xticklabels([]);

% Plot colored boxes around the significantly modulated timepoints
for neui = 1:nsig
    thissig = sigofem(neui,:);
    thishilo = sighilo(neui,:);
    starts = find(diff(thissig)==1);
    if thissig(1)==1
        starts = [0, starts];
    end
    stops = find(diff(thissig)==-1);
    if thissig(end)==1
        stops = [stops, length(thissig)];
    end
    for boxi = 1:length(starts)
        ishilo = thishilo(starts(boxi)+1);
        if ishilo == 1
            lstyle = '-r';
        elseif ishilo == -1
            lstyle = '-c';
        end
        lwidth = 1;
        plot([starts(boxi)+0.5 starts(boxi)+0.5], [neui-0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
        plot([stops(boxi)+0.5 stops(boxi)+0.5], [neui-0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
        plot([starts(boxi)+0.5 stops(boxi)+0.5], [neui-0.5 neui-0.5], lstyle, 'LineWidth', lwidth);
        plot([starts(boxi)+0.5 stops(boxi)+0.5], [neui+0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
    end
end

set(gca,'TickLength',[0.035 0.01]);

% Speed
figure; clf; hold on;
patch([plotx, fliplr(plotx)], [errdown, fliplr(errup)], 'k', 'FaceAlpha', 0.25, 'LineStyle', 'none');
plot(meanspeed, '-k', 'LineWidth', 4);
xticks([]);
set(gca, 'TickDir', 'out');
ax = gca;
ax.FontSize = 30;
xlim([0.5, npt+0.5]);
ylim([0 500]);
yticks([0 250 500]);
yticklabels([]);

plot([npt-65 npt-15], [150 150], '-k', 'LineWidth', 8)
set(gca,'TickLength',[0.03 0.01]);
set(gca, 'YAxisLocation', 'right');

% Sigs
figure; clf; hold on;
ntot = length(sigi);
plot(sum(sigofem,1)./ntot, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'black');
plot(sum(sighilo==-1,1)./ntot, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'cyan');
plot(sum(sighilo==1,1)./ntot, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'red');

xticks([0.5, Settings.Stretch.nbef+0.5-Settings.Stretch.periwin, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5+Settings.Stretch.periwin, npt+0.5]);
set(gca, 'XTickLabels', []);

a = get(gca, 'XTickLabels');
set(gca, 'XTickLabels', a, 'FontSize', 30);
set(gca, 'TickDir', 'out');

xlim([0.5, npt+0.5]);
set(gca,'TickLength',[0.03 0.01]);
yticklabels([]);
set(gca,'YAxisLocation', 'right');

% Find percentiles of modulation times
ntot = length(sigi);
timevec = -.395:.01:1.195;
sigsum = squeeze(sum(sigofem,1));
hisum = squeeze(sum(sighilo==1,1));
losum = squeeze(sum(sighilo==-1,1));

sigsum = sigsum(periwin);
hisum = hisum(periwin);
losum = losum(periwin);

sigcdf = cumsum(sigsum./ntot);
hicdf = cumsum(hisum./ntot);
locdf = cumsum(losum./ntot);

sigcdf = sigcdf./sigcdf(end);
hicdf = hicdf./hicdf(end);
locdf = locdf./locdf(end);

sigtimes = [];
hitimes = [];
lotimes = [];
for si = 1:length(sigsum)
    sigtimes = [sigtimes; ones(sigsum(si),1)*timevec(si)];
    hitimes = [hitimes; ones(hisum(si),1)*timevec(si)];
    lotimes = [lotimes; ones(losum(si),1)*timevec(si)];
end
medsig = median(sigtimes);
medhi = median(hitimes);
medlo = median(lotimes);

figure; clf; hold on;
axis square

startpctilesig = sum(sigtimes < 0)/numel(sigtimes);
endpctilesig = sum(sigtimes < 0.8)/numel(sigtimes);

plot([timevec(1) timevec(end)], [0 1], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 2);
plot([0.4 0.4], [0 1], '--', 'Color', [0.75 0.75 0.75]);
plot([timevec(1) timevec(end)], [0.5 0.5], '--', 'Color', [0.75 0.75 0.75]);

stairs(timevec, sigcdf, '-k', 'LineWidth', 2);
stairs(timevec, hicdf, '-r', 'LineWidth', 1);
stairs(timevec, locdf, '-c', 'LineWidth', 1);

% Plot out the percentiles but just for allsig
plot([medsig medsig], [0 0.5], '-', 'Color', 'k');
plot([0 0], [0 startpctilesig], '-', 'Color', 'k');
plot([-0.4 0], [startpctilesig startpctilesig], '-', 'Color', 'k');
plot([0.8 0.8], [0 endpctilesig], '-', 'Color', 'k');
plot([-0.4 0.8], [endpctilesig endpctilesig], '-', 'Color', 'k');
xlim([-.4 1.2]);
xticks([-.4 0 0.8 1.2]);
xticklabels([]);
yticks([0 startpctilesig 0.5 endpctilesig 1]);
yticklabels([]);
set(gca, 'XTickLabels', []);
a = get(gca, 'XTickLabels');
set(gca, 'XTickLabels', a, 'FontSize', 30);
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.03 0.01]);

R(alphi).startpctilesig = startpctilesig;
R(alphi).endpctilesig = endpctilesig;
R(alphi).medsig = medsig;
R(alphi).medhi = medhi;
R(alphi).medlo = medlo;
R(alphi).psigtimes = signrank(sigtimes-0.4);
R(alphi).phitimes = signrank(hitimes-0.4);
R(alphi).plotimes = signrank(lotimes-0.4);

%% now do locations
% CHANGE THIS FOR RELATIVE OR ABSOLUTE DEPTH
usedepths = Settings.Global.Locs.depths;

locs.center = 2;
locs.anterior = 1;
locs.posterior = 3;
locs.medial = 4;
orderednames = {'anterior', 'center', 'posterior', 'medial'};

lwid = 0.2;
lgap = 0.05;

% Now get where the sigs are and where the nonsigs are
for typei = 1:Settings.Global.Locs.ntracktypes
    tracksigdepths{typei} = [];
    tracknotsigdepths{typei} = [];
end

allsigdepths = Settings.Global.Locs.depths(sigi);
allnotsigdepths = Settings.Global.Locs.depths(~sigi);

for alli = 1:length(Settings.Global.Locs.depths)
    xloc = locs.(Settings.Global.Locs.tracks{alli});
    if sigi(alli)
        tracksigdepths{xloc} = [tracksigdepths{xloc}; Settings.Global.Locs.depths(alli)];
    else
        tracknotsigdepths{xloc} = [tracknotsigdepths{xloc}; Settings.Global.Locs.depths(alli)];
    end
end

% do stats
trackdepthpvals = nan(Settings.Global.Locs.ntracktypes,1);
for typei = 1:Settings.Global.Locs.ntracktypes
    if ~(isempty(tracksigdepths{typei}) || isempty(tracknotsigdepths{typei}))
        trackdepthpvals(typei) = ranksum(tracksigdepths{typei}, tracknotsigdepths{typei});
    end
end
depthpvals = ranksum(allsigdepths,allnotsigdepths);

% plot them

figure;
hold on;
for typei = 1:Settings.Global.Locs.ntracktypes
    thisloc = typei;
    plot(thisloc*[1 1], [0 13.4], '-k');
    plot(thisloc*ones(length(tracksigdepths{typei}),1)-lwid/2, tracksigdepths{typei}, '.r')
    plot([thisloc-lwid thisloc], median(tracksigdepths{typei})*[1 1], '-r', 'LineWidth', 3)
    plot(thisloc*ones(length(tracknotsigdepths{typei}),1)+lwid/2, tracknotsigdepths{typei}, '.k')
    plot([thisloc thisloc+lwid], median(tracknotsigdepths{typei})*[1 1], '-k', 'LineWidth', 3)
    if trackdepthpvals(typei) < 0.05
        text(thisloc, 14, ['p = ' num2str(trackdepthpvals(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'r');
    else
        text(thisloc, 14, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'k');
    end
end

thisloc = Settings.Global.Locs.ntracktypes+2;
plot(thisloc*ones(length(allsigdepths),1)-lwid/2, allsigdepths, '.r')
plot([thisloc-lwid thisloc], median(allsigdepths)*[1 1], '-r', 'LineWidth', 3)
plot(thisloc*ones(length(allnotsigdepths),1)+lwid/2, allnotsigdepths, '.k')
plot([thisloc thisloc+lwid], median(allnotsigdepths)*[1 1], '-k', 'LineWidth', 3)
if depthpvals < 0.05
    text(thisloc, 14, ['p = ' num2str(depthpvals)], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'r');
else
    text(thisloc, 14, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'k');
end

xticks([1:Settings.Global.Locs.ntracktypes Settings.Global.Locs.ntracktypes+2]);
xticklabels({'ant', 'ctr', 'post', 'med', 'ALL'});
grid on;
xlim([0 Settings.Global.Locs.ntracktypes+3]);
yticks([0 5 10]);
if all(usedepths == Settings.Global.Locs.depths)
    ylabel('Depth (absolute, mm)');
else
    ylabel('Depth (relative, mm)');
end

% what about track proportions
totprop = sum(sigi)/length(sigi);
trackprops = nan(Settings.Global.Locs.ntracktypes,1);
trackns = nan(Settings.Global.Locs.ntracktypes,1);
for typei = 1:Settings.Global.Locs.ntracktypes
    trackns(typei) = length(tracksigdepths{typei}) + length(tracknotsigdepths{typei});
    trackprops(typei) = length(tracksigdepths{typei})/trackns(typei);
end

npair = Settings.Global.Locs.ntracktypes*(Settings.Global.Locs.ntracktypes-1)/2;
pairsigsz = [];
pairsigschi = [];
for typei1 = 1:Settings.Global.Locs.ntracktypes-1
    for typei2 = typei1+1:Settings.Global.Locs.ntracktypes
        jointprop = (trackns(typei1)*trackprops(typei1) + trackns(typei2)*trackprops(typei2))/(trackns(typei1)+trackns(typei2));

        % Do it with z test
        tspz = (trackprops(typei1)-trackprops(typei2))/sqrt(jointprop*(1-jointprop)*((1/trackns(typei1))+(1/trackns(typei2))));
        if tspz > 0
            tspp = 2*(1-normcdf(tspz));
        else
            tspp = 2*normcdf(tspz);
        end
        if tspp < 0.05/npair
            pairsigsz = [pairsigsz; [typei1 typei2]];
        end

        % Do it with chi-sq test
        n1 = length(tracksigdepths{typei1});
        n2 = length(tracksigdepths{typei2});
        N1 = trackns(typei1);
        N2 = trackns(typei2);
        x1 = [repmat('a',N1,1); repmat('b',N2,1)];
        x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
        [tbl,chi2stat,chipval] = crosstab(x1,x2);
        if chipval < 0.05/npair
            pairsigschi = [pairsigschi; [typei1 typei2]];
        end
    end
end

ylabel([]);
yticklabels([]);
xticklabels([]);
title('');

figure('Renderer','Painters');
hold on;
for typei = 1:Settings.Global.Locs.ntracktypes
    bar(typei, trackprops(typei), 'EdgeColor', 'none', 'FaceColor', 'k')
    text(typei, 0, ['n=' num2str(trackns(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'w');
end
plot([0 Settings.Global.Locs.ntracktypes+1], [totprop totprop], '--r');
xlim([0 Settings.Global.Locs.ntracktypes+1]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xticks(1:Settings.Global.Locs.ntracktypes);
xticklabels({'ant', 'ctr', 'post', 'med'});
text(Settings.Global.Locs.ntracktypes+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'r');
if isempty(pairsigschi)
    text(mean([1 Settings.Global.Locs.ntracktypes]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
else
    for pairi = 1:size(pairsigschi,1)
        plot([pairsigschi(1) pairsigschi(2)], [1 1]-pairi*0.05, '-k');
    end
end

ylabel([]);
yticklabels([]);
xticklabels([]);
title('');

varargout{1} = R;
