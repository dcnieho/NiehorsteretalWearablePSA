% Output figures:
%   1) apparent gaze shift functions with shaded MAD regions
%   2) absolute slope dot plot by eye tracker for the center target only
%   3) signed slope plot: left eye / version / right eye for the center target
%   4) absolute slope bar graph: left / center / right targets
close all; clear; clc;

%% Configuration
dataFolder = fullfile('.', 'data');
plotFolder = fullfile('.', 'figures');
if ~exist(plotFolder, 'dir'); mkdir(plotFolder); end
refPupilSize = 4;
nBins = 11;
binWidth = 1/nBins;
binStart = binWidth/2;
binEnd = 1-binStart;
minDatToInclude = 2500; % ms per bin, as in the original analysis
maxMAD = 0.6;           % deg; reject participant/tracker/target if mean MAD exceeds this
windowLength = 50;      % ms, RIONEPS window
growWindow = 20;        % ms, expansion around RIONEPS detections

% Sampling rates, reference-offset sample counts, and RIONEPS thresholds.
fs = {'invisible',200; 'neon',200; 'smi',60; 'g2',50; 'g3',50};
nSamp = {'invisible',140; 'neon',140; 'smi',42; 'g2',35; 'g3',35};  % 700 ms
IMThresholds = {'invisible',[15 15 15]; 'neon',[300 300 300]; ...
                'smi',[8 8 8]; 'g2',[8 8 8]; 'g3',[8 8 8]};

% Display labels and colors
etLbls = {'Tobii Glasses 2','Tobii Glasses 3','Pupil Invisible','Pupil Neon','SMI ETG 2'};
etLblsnl = {'Tobii\newlineGlasses 2','Tobii\newlineGlasses 3', 'Pupil\newlineInvisible','Pupil\newlineNeon','SMI\newlineETG 2'};
[etLbls, order] = sort(etLbls);
etLblsnl = etLblsnl(order);
colors = [230,25,75; 255,225,25; 60,180,75; 0,130,200; 82,14,128; 253,95,22]./255;

%% 1. File selection and metadata
files = dir(fullfile(dataFolder, '*'));
files = files(~[files.isdir]);
dataFiles = struct('folder',{},'name',{},'participant',{},'eyeTracker',{},'target',{},'trial',{});
for i = 1:numel(files)
    [~, stem] = fileparts(files(i).name);
    parts = split(stem, '_');
    if numel(parts) < 4; continue; end
    dataFiles(end+1).folder = files(i).folder; %#ok<SAGROW>
    dataFiles(end).name = files(i).name;
    dataFiles(end).participant = char(parts(1));
    dataFiles(end).eyeTracker = char(parts(2));
    dataFiles(end).target = char(parts(3));
    dataFiles(end).trial = char(parts(4));
end
participants = sort(unique({dataFiles.participant}));
eyeTrackers = unique({dataFiles.eyeTracker});
eyeTrackers = eyeTrackers(order);
targets = unique({dataFiles.target});
trials = unique({dataFiles.trial});
assert(numel(participants)==6, 'Expected packaged participants P1-P6.');
assert(all(strcmp(participants, arrayfun(@(x)sprintf('P%d',x),1:6,'UniformOutput',false))), 'Packaged participant labels must be P1-P6.');
centerTarget = find(strcmp(targets,'center'),1);
leftTarget = find(strcmp(targets,'left'),1);
rightTarget = find(strcmp(targets,'right'),1);
assert(~isempty(centerTarget) && ~isempty(leftTarget) && ~isempty(rightTarget), 'Expected targets named left, center, and right.');

%% 2. Cleanup and data selection
% apparentGazeShift{tracker,participant,target,trial,signal}
% signal: 1 version-H, 2 version-V, 3 left-H, 4 left-V, 5 right-H, 6 right-V
apparentGazeShift = cell(numel(eyeTrackers),numel(participants),numel(targets),numel(trials),6);
for e = 1:numel(eyeTrackers)
    trackerMask = strcmp({dataFiles.eyeTracker}, eyeTrackers{e});
    tfs = lookupValue(fs, eyeTrackers{e});
    nOffsetSamples = lookupValue(nSamp, eyeTrackers{e});
    imThreshold = lookupValue(IMThresholds, eyeTrackers{e});
    for s = 1:numel(participants)
        participantMask = trackerMask & strcmp({dataFiles.participant}, participants{s});
        [blinkOn, blinkOff] = deal({[],[],[]});
        [saccadeOn, saccadeOff] = deal({[],[],[]});
        for t = 1:numel(targets)
            targetMask = participantMask & strcmp({dataFiles.target}, targets{t});
            for r = 1:numel(trials)
                fileMask = targetMask & strcmp({dataFiles.trial}, trials{r});
                if ~any(fileMask); continue; end
                assert(sum(fileMask)==1, 'Expected one file per tracker/participant/target/trial.');
                data = readtable(fullfile(dataFiles(fileMask).folder, dataFiles(fileMask).name), ...
                    'FileType','delimitedtext');
                meanPupil = mean([data.Basler_left_diameter data.Basler_right_diameter],2);

                % Reference each gaze signal to samples nearest 4-mm pupil size.
                [~, idx] = sort(abs(meanPupil-refPupilSize));
                idx = idx(1:min(nOffsetSamples,numel(idx)));
                offset = mean(data.offset_x_binoc(idx),'omitnan');
                apparentGazeShift{e,s,t,r,1} = [meanPupil data.offset_x_binoc-offset data.timestamp];
                offset = mean(data.offset_y_binoc(idx),'omitnan');
                apparentGazeShift{e,s,t,r,2} = [meanPupil -(data.offset_y_binoc-offset) data.timestamp];
                haveIndividualEyes = any(strcmp(data.Properties.VariableNames,'offset_x_left'));
                if haveIndividualEyes
                    [~, idx] = sort(abs(data.Basler_left_diameter-refPupilSize));
                    idx = idx(1:min(nOffsetSamples,numel(idx)));
                    offset = mean(data.offset_x_left(idx),'omitnan');
                    apparentGazeShift{e,s,t,r,3} = [data.Basler_left_diameter data.offset_x_left-offset];
                    offset = mean(data.offset_y_left(idx),'omitnan');
                    apparentGazeShift{e,s,t,r,4} = [data.Basler_left_diameter -(data.offset_y_left-offset)];
                    [~, idx] = sort(abs(data.Basler_right_diameter-refPupilSize));
                    idx = idx(1:min(nOffsetSamples,numel(idx)));
                    offset = mean(data.offset_x_right(idx),'omitnan');
                    apparentGazeShift{e,s,t,r,5} = [data.Basler_right_diameter data.offset_x_right-offset];
                    offset = mean(data.offset_y_right(idx),'omitnan');
                    apparentGazeShift{e,s,t,r,6} = [data.Basler_right_diameter -(data.offset_y_right-offset)];
                end

                firstTrial = find(~cellfun(@isempty, apparentGazeShift(e,s,t,:,1)),1);
                t0 = apparentGazeShift{e,s,t,firstTrial,1}(1,end);
                time = apparentGazeShift{e,s,t,r,1}(:,end)-t0;

                % ---- Blink removal ----
                filterWindow = 30; % ms
                window = ceil(filterWindow/1000*tfs);
                ntaps = 2*ceil(window)-1;
                pupilSize = apparentGazeShift{e,s,t,r,1}(:,1);
                pupilVelocity = -savitzkyGolayFilt(pupilSize,2,1,ntaps)*tfs;
                minFixSamples = ceil(50/1000*tfs);
                centralFixSamples = ceil(30/1000*tfs);
                exitCriterion = 5;
                nStd = 8;
                peakThreshold = 10;
                previousThreshold = inf;
                while previousThreshold-peakThreshold > exitCriterion
                    previousThreshold = peakThreshold;
                    qBelow = pupilVelocity < peakThreshold;
                    [on, off] = bool2bounds(qBelow);
                    qLong = off-on >= minFixSamples;
                    on = on(qLong) + centralFixSamples;
                    off = off(qLong) - centralFixSamples;
                    valid = off>=on;
                    on = on(valid); off = off(valid);
                    if isempty(on); break; end
                    idxFix = bounds2ind(on,off);
                    mu = mean(pupilVelocity(idxFix),'omitnan');
                    sigma = std(pupilVelocity(idxFix),'omitnan');
                    peakThreshold = mu+nStd*sigma;
                end
                peakThreshold = min(10,peakThreshold);
                qBlink = isnan(pupilVelocity) | pupilSize<1 | abs(pupilVelocity)>peakThreshold;
                [bon,boff] = bool2bounds(qBlink);
                dur = boff-bon+1;
                qRemoveBlink = dur<3 | bon==1 | boff==numel(pupilSize);
                bon(qRemoveBlink)=[]; boff(qRemoveBlink)=[];
                bon = max(1,bon-ceil(80/1000*tfs));
                boff = min(numel(pupilSize),boff+ceil(120/1000*tfs));
                blinkOn{t} = [blinkOn{t} time(bon).'];
                blinkOff{t} = [blinkOff{t} time(boff).'];
                qRemove = false(numel(pupilSize),1);
                % mark for removal
                for b = 1:numel(blinkOn{t})
                    qRemove = qRemove | (time>=blinkOn{t}(b) & time<=blinkOff{t}(b));
                end
                % and actually remove
                for sig = 1:(2+4*haveIndividualEyes)
                    apparentGazeShift{e,s,t,r,sig}(qRemove,1:2) = nan;
                end

                % ---- find and remove looks away from the target ----
                ek = EK([apparentGazeShift{e,s,t,r,1}(:,2) apparentGazeShift{e,s,t,r,2}(:,2)], ...
                    tfs,'','','','',5);
                if ~isempty(ek.ms)
                    ms = ek.ms;
                    match = find_opposite_saccades(ms,'MaxGap',ceil(400/1000*tfs));
                    if ~isempty(match)
                        sOn = match(:,3).'; sOff = match(:,6).';
                        ms([match(:,1);match(:,2)],:) = [];
                    else
                        sOn = []; sOff = [];
                    end
                    sOn = [sOn ms(:,1).']; sOff = [sOff ms(:,2).'];
                    saccadeOn{t} = [saccadeOn{t} time(sOn).'];
                    saccadeOff{t} = [saccadeOff{t} time(sOff).'];
                end
                for b = 1:numel(saccadeOn{t})
                    qRemove = qRemove | (time>=saccadeOn{t}(b) & time<=saccadeOff{t}(b));
                end
                % ---- remove spikes in the signal using RIONEPS detector ----
                qRemoveBySignal = {qRemove};
                if haveIndividualEyes; qRemoveBySignal = repmat(qRemoveBySignal,1,3); end
                growSamples = ceil(growWindow/1000*tfs);
                winSamples = ceil(windowLength/1000*tfs);
                ns = numel(time);
                for signalGroup = 1:(1+2*haveIndividualEyes)
                    pathLen = nan(ns,1); crowFliesDist = nan(ns,1);
                    for ii = 1:ns
                        if time(ii)-windowLength/1000/2 < time(1) || ...
                           time(ii)+windowLength/1000/2 > time(end)
                            continue
                        end
                        qDat = time>=time(ii)-windowLength/1000/2 & time<=time(ii)+windowLength/1000/2;
                        x = apparentGazeShift{e,s,t,r,(signalGroup-1)*2+1}(qDat,2);
                        y = apparentGazeShift{e,s,t,r,(signalGroup-1)*2+2}(qDat,2);
                        dx = diff(x); dy = diff(y);
                        pathLen(ii) = sum(hypot(dx,dy),'omitnan');
                        crowFliesDist(ii) = hypot(sum(dx,'omitnan'),sum(dy,'omitnan'));
                    end
                    inefficiency = (pathLen-crowFliesDist)*winSamples;
                    [ron,roff] = bool2bounds(inefficiency>imThreshold(signalGroup));
                    ron = max(ron-growSamples,1); roff = min(roff+growSamples,ns);
                    for m = numel(ron):-1:2
                        if ron(m)<roff(m-1)
                            roff(m-1)=roff(m); ron(m)=[]; roff(m)=[];
                        end
                    end
                    qRemoveBySignal{signalGroup} = qRemoveBySignal{signalGroup} | bounds2bool(ron,roff,ns).';
                end
                % Remove gaze samples farther than 7.5 deg from fixation.
                for signalGroup = 1:(1+2*haveIndividualEyes)
                    x = apparentGazeShift{e,s,t,r,(signalGroup-1)*2+1}(:,2);
                    y = apparentGazeShift{e,s,t,r,(signalGroup-1)*2+2}(:,2);
                    qRemoveBySignal{signalGroup} = qRemoveBySignal{signalGroup} | hypot(x,y)>7.5;
                end
                for sig = 1:(2+4*haveIndividualEyes)
                    apparentGazeShift{e,s,t,r,sig}(qRemoveBySignal{ceil(sig/2)},1:2) = nan;
                end
            end
        end
    end
end

%% 3. Binned apparent gaze shift functions and straight-line fits
nE = numel(eyeTrackers);
nS = numel(participants);
nT = numel(targets);

agsVersion = cell(nE,nS,nT,2);      % each cell: [pupil, shift, MAD, n]
agsEyes = cell(nE,nS,2,2);          % tracker, participant, eye(L/R), dimension(H/V), center target
slopesVersion = nan(nE,nS,nT,2);
slopesEyes = nan(nE,nS,2,2);

for e = 1:nE
    tfs = lookupValue(fs,eyeTrackers{e});
    minSamples = ceil(minDatToInclude/1000*tfs);

    for s = 1:nS
        % ---- Version/binocular signal: all targets ----
        for t = 1:nT
            plotDat = cell(1,2);
            plotX = cell(1,2);
            plotY = cell(1,2);
            plotE = cell(1,2);
            n = cell(1,2);
            qUse = cell(1,2);

            for dim = 1:2
                plotDat{dim} = cat(1,apparentGazeShift{e,s,t,:,dim});
                if isempty(plotDat{dim})
                    continue
                end

                [plotX{dim},plotY{dim},plotE{dim},n{dim},qUse{dim}] = ...
                    quantile2D(plotDat{dim}(:,1),plotDat{dim}(:,2), ...
                    binStart,binWidth,binEnd);
            end

            if isempty(plotDat{1}) || isempty(plotDat{2})
                continue
            end

            if mean(n{1}) < minSamples || mean(n{2}) < minSamples
                continue
            end

            if mean(plotE{1}) > maxMAD || mean(plotE{2}) > maxMAD
                continue
            end

            for dim = 1:2
                agsVersion{e,s,t,dim} = ...
                    [plotX{dim}(:) plotY{dim}(:) plotE{dim}(:) n{dim}(:)];
                qUseThis = and(qUse{:}) & ...
                    ~isnan(plotDat{dim}(:,1)) & ...
                    ~isnan(plotDat{dim}(:,2));

                if sum(qUseThis) >= 2
                    p = polyfit(plotDat{dim}(qUseThis,1), ...
                                plotDat{dim}(qUseThis,2),1);
                    slopesVersion(e,s,t,dim) = p(1);
                end
            end
        end

        % ---- Separate left/right-eye signals: center target only ----
        for eye = 1:2
            sig0 = 3 + (eye-1)*2;

            plotDat = cell(1,2);
            plotX = cell(1,2);
            plotY = cell(1,2);
            plotE = cell(1,2);
            n = cell(1,2);
            qUse = cell(1,2);

            for dim = 1:2
                sig = sig0 + dim - 1;
                plotDat{dim} = cat(1,apparentGazeShift{e,s,centerTarget,:,sig});

                if isempty(plotDat{dim})
                    continue
                end

                [plotX{dim},plotY{dim},plotE{dim},n{dim},qUse{dim}] = ...
                    quantile2D(plotDat{dim}(:,1),plotDat{dim}(:,2), ...
                    binStart,binWidth,binEnd);
            end

            if isempty(plotDat{1}) || isempty(plotDat{2})
                continue
            end

            if mean(n{1}) < minSamples || mean(n{2}) < minSamples
                continue
            end

            if mean(plotE{1}) > maxMAD || mean(plotE{2}) > maxMAD
                continue
            end

            for dim = 1:2
                agsEyes{e,s,eye,dim} = ...
                    [plotX{dim}(:) plotY{dim}(:) plotE{dim}(:) n{dim}(:)];

                qUseThis = and(qUse{:}) & ...
                    ~isnan(plotDat{dim}(:,1)) & ...
                    ~isnan(plotDat{dim}(:,2));

                if sum(qUseThis) >= 2
                    p = polyfit(plotDat{dim}(qUseThis,1), ...
                                plotDat{dim}(qUseThis,2),1);
                    slopesEyes(e,s,eye,dim) = p(1);
                end
            end
        end
    end
end

%% 4A. Apparent gaze shift functions with shaded MAD regions
xLimits = [2 6];
yLimits = [-3 3];
% ---- Center target, version signal ----
f = figure('Color','w');
tl = tiledlayout(2,nE,'TileIndexing','rowmajor','Padding','tight','TileSpacing','compact');
xlabel(tl,'Pupil size (mm)','FontSize',14,'FontWeight','bold');
ylabel(tl,'Apparent gaze shift (deg)','FontSize',14,'FontWeight','bold');
legendHandles = gobjects(1,nS);
for dim = 1:2
    for e = 1:nE
        ax = nexttile((dim-1)*nE+e); hold(ax,'on');
        yline(ax,0,'--','Color',[.6 .6 .6]);
        for s = nS:-1:1
            dat = agsVersion{e,s,centerTarget,dim};
            if isempty(dat); continue; end
            h = plotMadBand(ax,dat,colors(s,:));
            if dim==1 && e==1; legendHandles(s)=h; end
        end
        xlim(ax,xLimits); ylim(ax,yLimits);
        ax.XTick=2:6; ax.YTick=-3:3; ax.FontSize=10;
        if dim==1; title(ax,etLbls{e},'FontSize',11); end
        if e==1
            ylabel(ax,ternary(dim==1,'Horizontal','Vertical'),'FontWeight','bold');
        end
    end
end
q=isgraphics(legendHandles);
hl = legend(legendHandles(q),participants(q),'Location','southwest','NumColumns',2,'Box','off');
f.Position(3)=f.Position(3)*2.2; f.Position(4)=f.Position(4)*1.2;
exportgraphics(f,fullfile(plotFolder,'center_binoc_error.png'),'Resolution',300);
% ---- Left/right targets, version signal ----
% Outer rows: Horizontal / Vertical
% Inner rows: Left target / Right target
sideTargets = [leftTarget rightTarget];
sideLabels = {'Left target','Right target'};
f = figure('Color','w');
tl = tiledlayout(2,1, ...
    'TileIndexing','rowmajor', ...
    'Padding','tight', ...
    'TileSpacing','compact');

xlabel(tl,'Pupil size (mm)', ...
    'FontSize',16, ...
    'FontWeight','bold');

ylabel(tl,'Apparent gaze shift (deg)', ...
    'FontSize',16, ...
    'FontWeight','bold');
for dim = 1:2
    if dim == 1
        tlh = tiledlayout(tl,2,nE, ...
            'TileIndexing','rowmajor', ...
            'Padding','tight', ...
            'TileSpacing','compact');
        tlh.Layout.Tile = 1;

        ylabel(tlh,'Horizontal', ...
            'FontSize',14, ...
            'FontWeight','bold');
    else
        tlv = tiledlayout(tl,2,nE, ...
            'TileIndexing','rowmajor', ...
            'Padding','tight', ...
            'TileSpacing','compact');
        tlv.Layout.Tile = 2;

        ylabel(tlv,'Vertical', ...
            'FontSize',14, ...
            'FontWeight','bold');
    end

    for e = 1:nE
        for side = 1:2
            t = sideTargets(side);

            if dim == 1
                ax = nexttile(tlh,e + (side-1)*nE);
            else
                ax = nexttile(tlv,e + (side-1)*nE);
            end
            hold(ax,'on');

            yline(ax,0,'--','Color',[.6 .6 .6]);

            h = gobjects(1,nS);
            % Reverse order as in original so patches/lines overlay similarly
            for s = nS:-1:1
                dat = agsVersion{e,s,t,dim};

                if isempty(dat)
                    continue
                end

                h(s) = plotMadBand(ax,dat,colors(s,:));
            end

            % Legend comes from first panel
            if e==1 && dim==1 && side==1
                qUse = isgraphics(h);
                hl = legend(ax,h(qUse),participants(qUse), ...
                    'FontSize',10, ...
                    'Location','southwest', ...
                    'NumColumns',2, ...
                    'Box','off');
            end

            xlim(ax,xLimits);
            ylim(ax,yLimits);

            ax.XTick = 2:6;
            ax.YTick = -3:3;

            % Only bottom-most row gets x tick labels
            if ~(side==2 && dim==2)
                ax.XTickLabel = {};
            end

            % Target labels on left
            if e==1
                ylabel(ax,sideLabels{side}, ...
                    'FontSize',12, ...
                    'FontWeight','bold');
            end

            % Eye-tracker titles only across the top
            if side==1 && dim==1
                title(ax,etLbls{e},'FontSize',12);
            end

            ax.XAxis.FontSize = 11;
            ax.YAxis.FontSize = 11;
        end
    end
end
f.Position(3)=f.Position(3)*2; f.Position(4)=f.Position(4)*2.2; f.Position(1:2) = 100;
drawnow
hl.Position(2) = .52;
exportgraphics(f,fullfile(plotFolder,'left_right_binoc_v2_error.png'),'Resolution',300);

% ---- Center target, separate left/right-eye signals ----
% Outer rows: Horizontal / Vertical
% Inner rows: Left eye / Right eye
hasEyeData = squeeze(any(any(any(~isnan(slopesEyes),4),3),2));
eyeTrackerIdx = find(hasEyeData);
eyeLabels = {'Left eye','Right eye'};
nEyeTrackers = numel(eyeTrackerIdx);

f = figure('Color','w');
tl = tiledlayout(2,1, ...
    'TileIndexing','rowmajor', ...
    'Padding','tight', ...
    'TileSpacing','compact');

xlabel(tl,'Pupil size (mm)', ...
    'FontSize',16, ...
    'FontWeight','bold');

ylabel(tl,'Apparent gaze shift (deg)', ...
    'FontSize',16, ...
    'FontWeight','bold');
for dim = 1:2
    % Nested layout for horizontal / vertical
    if dim == 1
        tlh = tiledlayout(tl,2,nEyeTrackers, ...
            'TileIndexing','rowmajor', ...
            'Padding','tight', ...
            'TileSpacing','compact');
        tlh.Layout.Tile = 1;

        ylabel(tlh,'Horizontal', ...
            'FontSize',14, ...
            'FontWeight','bold');
    else
        tlv = tiledlayout(tl,2,nEyeTrackers, ...
            'TileIndexing','rowmajor', ...
            'Padding','tight', ...
            'TileSpacing','compact');
        tlv.Layout.Tile = 2;

        ylabel(tlv,'Vertical', ...
            'FontSize',14, ...
            'FontWeight','bold');
    end

    for ee = 1:nEyeTrackers
        e = eyeTrackerIdx(ee);

        for eye = 1:2
            if dim == 1
                ax = nexttile(tlh,ee + (eye-1)*nEyeTrackers);
            else
                ax = nexttile(tlv,ee + (eye-1)*nEyeTrackers);
            end
            hold(ax,'on');

            yline(ax,0,'--','Color',[.6 .6 .6]);

            h = gobjects(1,nS);

            for s = nS:-1:1
                dat = agsEyes{e,s,eye,dim};
                if isempty(dat)
                    continue
                end

                h(s) = plotMadBand(ax,dat,colors(s,:));
            end

            % Build legend from the upper-left panel.
            % Add off-axis dummy handles for participants without valid
            % data in this particular panel so P1-P6 always appear.
            if ee==1 && dim==1 && eye==1
                qMissing = ~isgraphics(h);

                xFake = xLimits(1) - diff(xLimits);
                yFake = yLimits(1) - diff(yLimits);

                for s = find(qMissing)
                    h(s) = plot(ax,xFake,yFake,'-', ...
                        'Color',colors(s,:), ...
                        'LineWidth',1.5);
                end

                hl = legend(ax,h,participants, ...
                    'FontSize',10, ...
                    'Location','southwest', ...
                    'NumColumns',2, ...
                    'Box','off');
            end

            xlim(ax,xLimits);
            ylim(ax,yLimits);

            ax.XTick = 2:6;
            ax.YTick = -3:3;

            % Only bottom-most row gets x tick labels
            if ~(eye==2 && dim==2)
                ax.XTickLabel = {};
            end

            % Eye labels on left
            if ee==1
                ylabel(ax,eyeLabels{eye}, ...
                    'FontSize',12, ...
                    'FontWeight','bold');
            end

            % Eye-tracker titles only across the top
            if eye==1 && dim==1
                title(ax,etLbls{e},'FontSize',12);
            end

            ax.XAxis.FontSize = 11;
            ax.YAxis.FontSize = 11;
        end
    end
end
f.Position(3)=f.Position(3)*1.8; f.Position(4)=f.Position(4)*2.2; f.Position(1:2) = 100;
drawnow
hl.Position(2) = .52;
exportgraphics(f,fullfile(plotFolder,'center_separate_eyes_error.png'),'Resolution',300);

%% 4B. Absolute slope dot plot by eye tracker, center target only
f = figure('Color','w');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
ylabel(tl,'Absolute slope (deg/mm)','FontSize',15,'FontWeight','bold');
participantOffset = linspace(-.11,.11,nS);
for dim = 1:2
    ax=nexttile; hold(ax,'on');
    h = gobjects(1,nS);
    dat = abs(squeeze(slopesVersion(:,:,centerTarget,dim))); % tracker x participant
    for s=1:nS
        h(s)=plot(ax,(1:nE)+participantOffset(s),dat(:,s),'o', ...
            'Color',colors(s,:),'MarkerFaceColor',colors(s,:),'MarkerSize',4);
    end
    mn=mean(dat,2,'omitnan');
    se=std(dat,0,2,'omitnan')./sqrt(sum(~isnan(dat),2));
    errorbar(ax,(1:nE)+.24,mn,se,'ko','MarkerFaceColor','k','LineWidth',1.5,'CapSize',0);
    title(ax,ternary(dim==1,'Horizontal','Vertical'),'FontWeight','bold');
    xlim(ax,[.7 nE+.3]); ylim(ax,[0 1.8]);
    ax.XTick=1:nE; ax.XTickLabel=etLblsnl; ax.FontSize=11;
    xlabel(ax,'Eye tracker','FontWeight','bold');
    if dim==1; legend(h,participants,'Location','northwest','NumColumns',2,'Box','off'); end
end
f.Position(3)=f.Position(3)*2.1;
exportgraphics(f,fullfile(plotFolder,'absolute_slope_dot_center.png'),'Resolution',300);

%% 4C. Signed slopes: left eye, version, right eye (center target)
hasEyeData = squeeze(any(any(any(~isnan(slopesEyes),4),3),2));
eyeTrackerIdx = find(hasEyeData);
f=figure('Color','w');
tl=tiledlayout(2,numel(eyeTrackerIdx),'TileIndexing','rowmajor','TileSpacing','compact','Padding','tight');
ylabel(tl,'Slope (deg/mm)','FontSize',15,'FontWeight','bold');
xlabel(tl,'Signal','FontSize',15,'FontWeight','bold');
signalLabels={'Left eye','Version','Right eye'}; xSignal=1:3;
legendHandles=gobjects(1,nS);
legendAx = gobjects(1);
for dim=1:2
    for ee=1:numel(eyeTrackerIdx)
        e=eyeTrackerIdx(ee); ax=nexttile((dim-1)*numel(eyeTrackerIdx)+ee); hold(ax,'on');
        if dim==1 && ee==1
            legendAx = ax;
        end
        yline(ax,0,'--','Color',[.6 .6 .6]);
        for s=1:nS
            dat=[slopesEyes(e,s,1,dim), slopesVersion(e,s,centerTarget,dim), slopesEyes(e,s,2,dim)];
            if any(isnan(dat)); continue; end
            plot(ax,xSignal,dat,'-o','Color',colors(s,:),'MarkerFaceColor',colors(s,:),'MarkerSize',6);
        end
        xlim(ax,[.5 3.5]); ylim(ax,[-2 1.2]); ax.XTick=xSignal; ax.XTickLabel=signalLabels;
        if dim==1; title(ax,etLbls{e},'FontWeight','bold'); ax.XTickLabel={}; end
        if ee==1; ylabel(ax,ternary(dim==1,'Horizontal','Vertical'),'FontWeight','bold'); end
        ax.FontSize=10;
    end
end
% Make legend entries independently of data availability so that
% all six participants P1-P6 are always represented.
for s = 1:nS
    legendHandles(s) = plot(legendAx,nan,nan,'-o', ...
        'Color',colors(s,:), ...
        'MarkerFaceColor',colors(s,:), ...
        'MarkerSize',6, ...
        'LineWidth',1);
end
legend(legendHandles,participants,'Location','southwest','NumColumns',2,'Box','off');
f.Position(3)=f.Position(3)*2.2; f.Position(4)=f.Position(4)*1.5;
exportgraphics(f,fullfile(plotFolder,'signed_slope_version_left_right_eye_center.png'),'Resolution',300);

%% 4D. Absolute slope bar graph: left, center, right targets
plotTargetIdx=[leftTarget centerTarget rightTarget];
plotTargetLabels={'Left','Center','Right'};
f=figure('Color','w');
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
ylabel(tl,'Absolute slope (deg/mm)','FontSize',15,'FontWeight','bold');
for dim=1:2
    ax=nexttile; hold(ax,'on');
    dat=abs(slopesVersion(:,:,plotTargetIdx,dim)); % tracker x participant x target
    mn=squeeze(mean(dat,2,'omitnan'));
    se=squeeze(std(dat,0,2,'omitnan'))./squeeze(sqrt(sum(~isnan(dat),2)));
    b=bar(ax,mn,'grouped');
    for k=1:3
        errorbar(ax,b(k).XEndPoints,mn(:,k),se(:,k),'k','LineStyle','none','LineWidth',1,'CapSize',8);
    end
    title(ax,ternary(dim==1,'Horizontal','Vertical'),'FontWeight','bold');
    xlim(ax,[.5 nE+.5]); ylim(ax,[0 1.2]);
    ax.XTick=1:nE; ax.XTickLabel=etLblsnl; ax.FontSize=11;
    xlabel(ax,'Eye tracker','FontWeight','bold');
    if dim==1; legend(ax,plotTargetLabels,'Location','northwest','Box','off'); end
end
f.Position(3)=f.Position(3)*2.1;
exportgraphics(f,fullfile(plotFolder,'absolute_slope_left_center_right.png'),'Resolution',300);

%% Local helper functions
function value = lookupValue(tableLike,key)
    idx=strcmp(tableLike(:,1),key);
    assert(any(idx),'No configuration found for eye tracker "%s".',key);
    value=tableLike{idx,2};
end

function h = plotMadBand(ax,dat,color)
    % dat columns: pupil size, median apparent gaze shift, MAD, n
    x = dat(:,1);
    y = dat(:,2);
    e = dat(:,3);
    good = ~isnan(x) & ~isnan(y) & ~isnan(e);
    x=x(good); y=y(good); e=e(good);
    if isempty(x)
        h = gobjects(1);
        return
    end
    xPatch = [x(:); flipud(x(:))];
    yPatch = [y(:)-e(:); flipud(y(:)+e(:))];
    patch(ax,xPatch,yPatch,color, ...
        'FaceAlpha',0.25, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');
    h = plot(ax,x,y,'-', ...
        'Color',color, ...
        'LineWidth',1.5);
end
function out=ternary(condition,a,b)
    if condition; out=a; else; out=b; end
end
