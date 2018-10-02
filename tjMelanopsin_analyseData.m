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
basePath = 'data';

%% Define the trial order
% Define the orders
theOrders = [3, 2, 1, 4 ; ...
    1, 4, 2, 3 ; ...
    1, 2, 3, 4 ; ...
    3, 1, 4, 2 ; ...
    1, 3, 2, 4 ; ...
    4, 2, 3, 1 ; ...
    1, 4, 3, 2 ; ...
    2, 3, 4, 1 ; ...
    3, 4, 2, 1 ; ...
    2, 4, 3, 1 ; ...
    4, 3, 2, 1 ; ...
    3, 1, 4, 2 ; ...
    2, 1, 4, 3 ; ...
    2, 4, 3, 1 ; ...
    4, 2, 3, 1 ; ...
    3, 4, 1, 2];

% Define the labels
theLabels = {'LFX' 'MEL' 'LMS' 'REF'};
theLabelsLong = {'Light flux' 'Melanopsin' 'LMS' 'Reference'};
theRGB = [0 0 0 ; ...
    67 162 202 ; ...
    254 196 79 ; ...
    189 189 189];

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
        tmpAnnot = dir(fullfile(theFolder, '0-*'));
        
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
        dataPath = fullfile(theFolder, '000', 'pupil_positions.csv');
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
        h.mainLine.Color = theRGB(ij, :)/255;
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
        w = [sin(2*pi*f*x) ; cos(2*pi*f*x)]' \ theQs;
        plot(x, [sin(2*pi*f*x) ; cos(2*pi*f*x)]'*w, '-r', 'LineWidth', 2); hold on;
        %plot(x, theQs);
        
        sw(ii, ij) = w(1);
        cw(ii, ij) = w(2);
        
        theAmplitude(ii, ij) = sw(ii, ij)*sqrt(1+((cw(ii, ij)/sw(ii, ij))^2));
        thePhase(ii, ij) = atan2(cw(ii, ij),sw(ii, ij));
        
        theComplexNum(ii, ij) = theAmplitude(ii, ij)*cos(thePhase(ii, ij)) + sqrt(-1)*theAmplitude(ii, ij)*sin(thePhase(ii, ij))
        c = c+1;
    end
    
end

set(sparkLineFigure, 'PaperPosition', [0 0 10 25]);
set(sparkLineFigure, 'PaperSize', [10 25]);
set(sparkLineFigure, 'Color', 'w');
set(sparkLineFigure, 'InvertHardcopy', 'off');
saveas(sparkLineFigure, 'Fig2.pdf', 'pdf');

%%
theAmplitudeMean = abs(mean(theComplexNum))
thePhaseMean = angle(mean(theComplexNum))



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

%%
polarPlotFigure = figure;
for ii = 1:4
    subplot(1, 5, ii);
    polar(0, 9, 'ok'); hold on;
    h0 = polar(thePhase(:, ii), 100*theAmplitude(:, ii), 'ok'); hold on
    h0.MarkerFaceColor = theRGB(ii, :)/255;
end
subplot(1, 5, 5);
for ii = 1:4
    h1 = polar(thePhaseMean(ii), 100*theAmplitudeMean(ii), 'ok'); hold on
    h1.MarkerSize = 10;
    h1.MarkerFaceColor = theRGB(ii, :)/255;
    h1.MarkerEdgeColor = 'r';
end

set(polarPlotFigure, 'PaperPosition', [0 0 30 10]);
set(polarPlotFigure, 'PaperSize', [30 10]);
set(polarPlotFigure, 'Color', 'w');
set(polarPlotFigure, 'InvertHardcopy', 'off');
saveas(polarPlotFigure, 'Fig3.pdf', 'pdf');

%% Plot
subplot(1, 2, 1);
plot([0 10], [0 10], ':k'); hold on
plot(100*abs(theComplexNum(:, 1)), 100*abs(theComplexNum(:, 3)+theComplexNum(:, 2)), 'ok', 'MarkerFaceColor', [44 162 95]/255);
xlim([0 10]); ylim([0 10]);
pbaspect([1 1 1]);
box off; set(gca, 'TickDir', 'out');
xlabel('Amplitude light flux [%]');
ylabel('Amplitude (LMS+Mel) [%]');
[RHO1,PVAL1] = corr(abs(theComplexNum(:, 1)), abs(theComplexNum(:, 4)))
[RHO2,PVAL2] = corr(abs(theComplexNum(:, 1)), abs(theComplexNum(:, 3)+theComplexNum(:, 2)))
text(4, 1, ['r=' num2str(RHO2, '%.2f') ', p=' num2str(PVAL2, '%.4f')]);
%text(4, 2, ['r=' num2str(RHO2) ', p=' num2str(PVAL2)]);

subplot(1, 2, 2);
plot([-pi pi], [-pi pi], ':k'); hold on
plot(angle(theComplexNum(:, 1)), angle(theComplexNum(:, 3)+theComplexNum(:, 2)), 'ok', 'MarkerFaceColor', [44 162 95]/255);
xlim([-pi pi]); ylim([-pi pi]);
pbaspect([1 1 1]);
box off; set(gca, 'TickDir', 'out');
xlabel('Phase light flux');
ylabel('Phase (LMS+Mel)');

set(gcf, 'PaperPosition', [0 0 15 10]);
set(gcf, 'PaperSize', [15 10]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Fig4.pdf', 'pdf');

%%
subplot(1, 3, 1);
plot(-theAmplitude(:, 1), -theAmplitude(:, 3), 'ok', 'MarkerFaceColor', 'k')
xlim([0 8]);
ylim([0 8]);
xlabel('Light flux Amplitude [%]');
ylabel('LMS Amplitude [%]');

subplot(1, 3, 2);
plot(-theAmplitude(:, 1), -theAmplitude(:, 2), 'ok', 'MarkerFaceColor', 'k')
xlim([0 8]);
ylim([0 8]);
xlabel('Light flux Amplitude [%]');
ylabel('Melanopsin Amplitude [%]');

subplot(1, 3, 3);
plot(-theAmplitude(:, 3), -theAmplitude(:, 2), 'ok', 'MarkerFaceColor', 'k')
xlim([0 8]);
ylim([0 8]);
xlabel('LMS Amplitude [%]');
ylabel('Melanopsin Amplitude [%]');


%% FFT
load('contrastTimeCourse.mat')

theAgeGroup = 1; % 1 = 20-24, 2 = 25-29, 3 = 30-35, 4 = 35-39

startIdx = 1202;

theRGB = SSTDefaultReceptorColors;

subplot(5, 4, 4);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).ref.lConeMean);
plot(f, xft, '-', 'Color', theRGB(1, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 1);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLFX.lConeMean);
plot(f, xft, '-', 'Color', theRGB(1, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 2);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modMel.lConeMean);
plot(f, xft, '-', 'Color', theRGB(1, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 3);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLMS.lConeMean);
plot(f, xft, '-', 'Color', theRGB(1, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 8);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).ref.mConeMean);
plot(f, xft, '-', 'Color', theRGB(2, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 5);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLFX.mConeMean);
plot(f, xft, '-', 'Color', theRGB(2, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 6);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modMel.mConeMean);
plot(f, xft, '-', 'Color', theRGB(2, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 7);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLMS.mConeMean);
plot(f, xft, '-', 'Color', theRGB(2, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 12);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).ref.sConeMean);
plot(f, xft, '-', 'Color', theRGB(3, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 9);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLFX.sConeMean);
plot(f, xft, '-', 'Color', theRGB(3, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 10);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modMel.sConeMean);
plot(f, xft, '-', 'Color', theRGB(3, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 11);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLMS.sConeMean);
plot(f, xft, '-', 'Color', theRGB(3, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 16);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).ref.melMean);
plot(f, xft, '-', 'Color', theRGB(4, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 13);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLFX.melMean);
plot(f, xft, '-', 'Color', theRGB(4, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 14);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modMel.melMean);
plot(f, xft, '-', 'Color', theRGB(4, :)); xlim([0.05 1]); ylim([0 1500]);
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

subplot(5, 4, 15);
[f, xft] = tjMelanopsin_doFFT(imgContent(theAgeGroup).modLMS.melMean);
plot(f, xft, '-', 'Color', theRGB(4, :)); xlim([0.05 1]); ylim([0 1500]); hold on;
xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
hold on; axis off; plot([0.05 0.05], [0 1500], '-k'); plot([0.05 1], [0 0], '-k');

cD = 1;
for ii = 1:4
    theData = [];
    for ij = 1:11
        theData = [theData theFFTPerCondition{ij, ii}];
    end
    subplot(5, 4, 16+ii);

    xlim([0.05 1]); ylim([0 500]);
    plot(0.5*[0.25 0.25], [0 500], ':r', 'LineWidth', 1.5);  hold on;
    plot([0.25 0.25], [0 500], ':r', 'LineWidth', 1.5);
    plot(2*[0.25 0.25], [0 500], ':r', 'LineWidth', 1.5);
        plot(fx, mean(theData, 2), '-k');
    xlabel('Frequency [Hz]'); pbaspect([1 1 1]); box off; ylabel('Power [aub]'); set (gca, 'TickDir', 'out');
    set(gca,'Color',[0.9 0.9 0.9])
    hold on; axis off; plot([0.05 0.05], [0 500], '-k'); plot([0.05 1], [0 0], '-k');
end
set(gcf, 'PaperPosition', [0 0 10 15]);
set(gcf, 'PaperSize', [10 15]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Fig5.pdf', 'pdf');