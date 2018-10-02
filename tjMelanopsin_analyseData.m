clearvars; close all;

%% Define the participant list
theSubjects = {'2018_02_06_s001' ...
    '2018_02_06_s002' ...
    '2018_02_06_s003' ...
    '2018_02_06_s004' ...
    '2018_02_07_s005' ...
    '2018_02_07_s006' ...
    '2018_02_07_s007' ...
    '2018_02_08_s008' ...
    '2018_02_09_s009' ...
    '2018_02_09_s010' ...
    '2018_02_09_s011'};

NSubjects = length(theSubjects);

% Define some paths
basePath = 'data';

%% Define the trial order
% Define the orders
theOrders = [1, 2, 3, 4 ; ...
    3, 4, 2, 1 ; ...
    3, 2, 1, 4 ; ...
    1, 3, 4, 2 ; ...
    3, 1, 2, 4 ; ...
    4, 2, 1, 3 ; ...
    3, 4, 1, 2 ; ...
    2, 1, 4, 3 ; ...
    1, 4, 2, 3 ; ...
    2, 4, 1, 3 ; ...
    4, 1, 2, 3 ; ...
    1, 3, 4, 2 ; ...
    2, 3, 4, 1 ; ...
    2, 4, 1, 3 ; ...
    4, 2, 1, 3 ; ...
    1, 4, 3, 2];



% Define the labels
theLabels = {'LMS' 'MEL' 'LFX' 'REF'}; 
theLabelsLong = {'LMS'  'Melanopsin' 'Light flux' 'Reference'};
theRGB = [254 196 79 ; ...
    67 162 202 ; ...
    10 10 10 ; ...
    255 255 255]/255;

%% Define camera info
frequencyHz = 30;
dt = 1/frequencyHz; % Camera sampling in Hz
endIdx = 30000;
x = 0:dt:(dt*(endIdx-1));
stimFreqHz = 0.25;

NTotalCycles = 90;
cycleStepSamples = 120;
theCyclesStart0 = 1:cycleStepSamples:(cycleStepSamples*NTotalCycles);
theCyclesEnd0 = cycleStepSamples:cycleStepSamples:(cycleStepSamples*NTotalCycles);

sparkLineFigure = figure;
c = 1;

%% Define participant
for ii = 1:NSubjects
    participantNum = ii;
    
    % Pick the order based on the input
    theOrderHere = theOrders(participantNum, :);
    
    for ij = 1:4
        theFolder = fullfile(basePath, theSubjects{participantNum}, num2str(ij-1, '%03g'), 'exports');
        
        %% Load the table
        tmpAnnot = dir(fullfile(theFolder, '0*'));
        
        %% Load the annotations
        % Annotation indices are:
        %   1 screen off
        %   2 start WB trailer
        %   3 end WB trailer
        %   4 start MGM trailer
        %   5 title slide
        %   6 end trailer (T&J)
        annotPath = fullfile(theFolder, tmpAnnot.name, 'annotations.csv');
        annotIdx = tjMelanopsin_loadAnnotationFile(annotPath);
        
        %% Load the data
        dataPath = fullfile(theFolder, tmpAnnot.name, 'pupil_positions.csv');
        [dataTraceRaw, timeTraceRaw, dataTraceIdx] = tjMelanopsin_loadDataFile(dataPath);
        
        % Associate the indices
        eyeMovementIdx = tjMelanopsin_associateIndices(dataTraceIdx, annotIdx);
        dataTraceRawExpt = dataTraceRaw(eyeMovementIdx(5):eyeMovementIdx(6));
        timeTraceRawExpt = timeTraceRaw(eyeMovementIdx(5):eyeMovementIdx(6));
        timeTraceRawExpt = timeTraceRawExpt-timeTraceRawExpt(1);
        
        % Process the data
        dataTraceRawExpt = tjMelanopsin_removeMissingData(dataTraceRawExpt);
        dataTraceInterpolated = inpaint_nans(dataTraceRawExpt);
        [dataTraceInterpolated, t] = tjMelanopsin_interpolateData(timeTraceRawExpt, dataTraceInterpolated, frequencyHz);
        dataTraceFiltered = tjMelanopsin_filterData(dataTraceInterpolated);
        dataTraceFiltered = dataTraceFiltered(1:9700);
        
        t = 0:(1/30):((length(dataTraceFiltered)-1)*(1/30));
        [pxx, fx] = plomb(dataTraceFiltered, t');
        %figure;
        %plot(f, xft); xlim([0.05 2]); ylim([0 600]);
        
        %% Average per cycle
        validCycles = length(dataTraceFiltered) > theCyclesEnd0;
        theCyclesStart = theCyclesStart0(validCycles);
        theCyclesEnd = theCyclesEnd0(validCycles);
        
        tmp = [];
        for im = 1:length(theCyclesStart)
            tmp = [tmp ; (dataTraceFiltered(theCyclesStart(im):theCyclesEnd(im))-nanmean(dataTraceFiltered))./nanmean(dataTraceFiltered)];
        end
        theDataPerCondition{theOrderHere(ij)} = tmp';
        theFFTPerCondition{ii, theOrderHere(ij)} = pxx;
    end
    
    %%
    figure(sparkLineFigure);
    dt = 1/30;
    x = 0:dt:(1/stimFreqHz-dt);
    for ij = 1:4
        theMeanAggregate{ii, ij} = mean(theDataPerCondition{ij}, 2);
        theSD = std(theDataPerCondition{ij}, [], 2);
        theSEM = theSD/sqrt(size(theDataPerCondition{ij}, 2));
        
        % Scale to be % change
        theQs = theMeanAggregate{ii, ij};
        %theQs = (theQ-mean(theQ))/mean(theQ);
        
        subplot_tight(11, 4, c);
        h = shadedErrorBar(x, theQs, theSEM); hold on;
        h.mainLine.Color = theRGB(ij, :);
        h.mainLine.LineWidth = 2;
        h.edge(1).LineStyle = 'none';
        h.edge(2).LineStyle = 'none';
        h.patch.FaceColor = [0.75 0.75 0.75];
        plot([x(1) x(end)], [0 0], '-k');
        plot([x(1) x(1)], [-0.1 0.1], '-k');
        
        %xlabel('Time [sec]');
        %ylabel('Pupil diameter [mm]');
        xlim([-0.1 4.1]);
        ylim([-0.1 0.1]);
        pbaspect([1 1 1]);
        box off;
        axis off;
        
        % Fit sine and cosine
        f = 0.25;
        w = [sin(2*pi*f*x) ; cos(2*pi*f*x)]' \ -theQs;
        plot(x, -[sin(2*pi*f*x) ; cos(2*pi*f*x)]'*w, '-r', 'LineWidth', 2); hold on;
        %plot(x, theQs);
        
        sw(ii, ij) = w(1);
        cw(ii, ij) = w(2);
        
        theAmplitude(ii, ij) = sqrt((cw(ii, ij).^2+sw(ii, ij).^2));
        thePhase(ii, ij) = atan2(cw(ii, ij),sw(ii, ij));
        
        theComplexNum(ii, ij) = theAmplitude(ii, ij)*cos(thePhase(ii, ij)) + sqrt(-1)*theAmplitude(ii, ij)*sin(thePhase(ii, ij))
        c = c+1;
    end
    
end

set(sparkLineFigure, 'PaperPosition', [0 0 10 25]);
set(sparkLineFigure, 'PaperSize', [10 25]);
set(sparkLineFigure, 'Color', 'w');
set(sparkLineFigure, 'InvertHardcopy', 'off');
saveas(sparkLineFigure, 'figures/Fig2.pdf', 'pdf');

%%
theAmplitudeMean = abs(mean(theComplexNum))
thePhaseMean = angle(mean(theComplexNum))

%% Show amplitude and phase
ampPhaseFigure = figure;
theAmplitudePct = theAmplitude*100;
subplot(2, 1, 1);
for ii = 1:4
plot((theAmplitudePct(:, ii)), 5-ii, 'ok', 'MarkerFaceColor', theRGB(ii, :)); hold on
theMean = mean(theAmplitudePct(:, ii));
theSD = std((theAmplitudePct(:, ii)));
theSEM = theSD / sqrt(size(theAmplitudePct, 1));
plot([theMean theMean], [5-(ii-0.3) 5-(ii+0.3)], '-r', 'LineWidth', 1.2); hold on
plot([theMean-theSEM theMean+theSEM], [5-ii 5-ii], '-r', 'LineWidth', 1.2); hold on
if ii == 4
   plot([theMean theMean], [0 5], ':k'); 
end
end
xlim([0 10]);
ylim([0 5]);
pbaspect([1 0.4 1]);
xlabel('Amplitude [\Delta%]');
ylabel('Direction');
set(gca, 'YTick', 1:4, 'YTickLabel', {theLabelsLong{end:-1:1}});
%set(gca, 'XTick', [-180:60:180]);
plot([0 0], [0 5], ':k');
box off; set(gca, 'TickDir', 'out');
title('Amplitude');


subplot(2, 1, 2);
for ii = 1:4
plot(rad2deg(thePhase(:, ii)), 5-ii, 'ok', 'MarkerFaceColor', theRGB(ii, :)); hold on
theMean = rad2deg(circ_mean((thePhase(:, ii))));
theSD = rad2deg(circ_std((thePhase(:, ii))));
theSEM = theSD / sqrt(size(thePhase, 1));
plot([theMean theMean], [5-(ii-0.3) 5-(ii+0.3)], '-r', 'LineWidth', 1.2); hold on
plot([theMean-theSEM theMean+theSEM], [5-ii 5-ii], '-r', 'LineWidth', 1.2); hold on
end
xlim([-180 180]); ylim([0 5]);
pbaspect([1 0.4 1]);
xlabel('Phase angle [\circ]');
ylabel('Direction');
set(gca, 'YTick', 1:4, 'YTickLabel', {theLabelsLong{end:-1:1}});
set(gca, 'XTick', [-180:60:180]);
plot([0 0], [0 5], ':k');
box off; set(gca, 'TickDir', 'out');
title('Phase');


set(ampPhaseFigure, 'PaperPosition', [0 0 10 10]);
set(ampPhaseFigure, 'PaperSize', [10 10]);
set(ampPhaseFigure, 'Color', 'w');
set(ampPhaseFigure, 'InvertHardcopy', 'off');
saveas(ampPhaseFigure, 'figures/Fig3.pdf', 'pdf');

%% Calculate the difference
for ii = 1:11
    subplot(11, 1, ii);
    plot([x(1) x(end)], [0 0], '-k'); hold on
    plot([x(1) x(1)], [-0.1 0.1], '-k');
    plot(x, theMeanAggregate{ii, 1} - theMeanAggregate{ii, 3}), '-r'; hold on
    plot(x, theMeanAggregate{ii, 2}, '-k');
    xlim([-0.1 4.1]);
    ylim([-0.1 0.1]);
    pbaspect([1 1 1]);
    box off;
    axis off;
end


%% Plot
summationFigure = figure;
plot([0 10], [0 10], ':k'); hold on
plot(100*abs(theComplexNum(:, 3)), 100*abs(theComplexNum(:, 1)+theComplexNum(:, 2)), 'ok', 'MarkerFaceColor', [44 162 95]/255);
plot(100*abs(theComplexNum(:, 3)), 100*abs(theComplexNum(:, 1)), 'sk', 'MarkerFaceColor', 0.8*[254 196 79]/255, 'MarkerEdgeColor', 'k');
xlim([0 10]); ylim([0 10]);
pbaspect([1 1 1]);
box off; set(gca, 'TickDir', 'out');
xlabel('Amplitude light flux [%]');
ylabel('Amplitude composite [%]');
disp('LMS and light flux')
[RHO1,PVAL1] = corr(abs(theComplexNum(:, 3)), abs(theComplexNum(:, 1)))

disp('LMS+mel and light flux')
[RHO2,PVAL2] = corr(abs(theComplexNum(:, 3)), abs(theComplexNum(:, 1)+theComplexNum(:, 2)))
text(4, 1, ['LMS+Mel: r=' num2str(RHO2, '%.2f') ', p=' num2str(PVAL2, '%.4f')]);
text(4, 2, ['LMS only: r=' num2str(RHO1, '%.2f') ', p=' num2str(PVAL1, '%.4f')]);


set(summationFigure, 'PaperPosition', [0 0 10 10]);
set(summationFigure, 'PaperSize', [10 10]);
set(summationFigure, 'Color', 'w');
set(summationFigure, 'InvertHardcopy', 'off');
saveas(summationFigure, 'figures/Fig4.pdf', 'pdf');
