function getStdOfRealData(roiToView)

cd ~/data/emri/paperdata
SIDs = ["s0615" "s0616" "s0617" "s0618" "s0619"];

allVoxelStds = [];
avgAllVoxelStdEachScan = [];
voxelAverageExponentImprovementEachScan = [];
voxelAverageRatioImprovementEachScan = [];
numAverages = [];
voxelStdOverNumAverages = [];


%get all of the data
for subject = 1:length(SIDs);

    cd(SIDs(subject));
    %create the view and go to unaveraged ungated scans
    v = newView;
    v = viewSet(v,'curGroup','rawUngatedResized')
    numScans = viewGet(v,'numScans');
    %define the roi time series that you want 
    %v = loadROI(v,'horizontalROI') ;
    %loop through each scan, getting individual and averaged std of tseries for the roi
    for scan = 1:numScans
        %set the current scan
        v = viewSet(v,'curScan',scan)    
        %load the t series
        roi = loadROITSeries(v,roiToView,viewGet(v,'curScan'),viewGet(v,'curGroup'));

        %get the std of all the voxels in the scan
            alltSeriesMeans = mean(roi.tSeries,2);
            alltSeries = 100*(roi.tSeries-alltSeriesMeans)./repmat(alltSeriesMeans,[1 size(roi.tSeries,2)]);
            allVoxelStds = [allVoxelStds std(alltSeries,0,2)'];
            avgSingleVoxelStd = mean(std(alltSeries,0,2));

        %get the std of all of the voxels averaged together
        tSeries = roi.tSeries;
        tSeries = mean(tSeries,1)';
        %convert to % signal change
        meanTSeries = mean(tSeries);
        tSeries = 100*(tSeries-meanTSeries)/meanTSeries;
        %save the std of the averaged voxels in this scan
        avgAllVoxelStdEachScan = [avgAllVoxelStdEachScan std(tSeries)];

        %check the improvement averaging over each scan
        [exponent, ratio, numAverages, voxelStdOverNumAverages] = checkAverageNoiseReduction(roi, avgSingleVoxelStd, subject, scan, numAverages, voxelStdOverNumAverages);
        voxelAverageExponentImprovementEachScan = [voxelAverageExponentImprovementEachScan exponent];
        voxelAverageRatioImprovementEachScan = [voxelAverageRatioImprovementEachScan ratio];

    end
    %go back and do it over
    pause(1)
    deleteView(v)
    mrQuit
    cd ..
   
end

%plot the improvement with averaging more voxels
figure(155), subplot(1,2,2)
hist(voxelAverageRatioImprovementEachScan); xlim([0 1])
xlabel('Ratio of observed to expected (1/sqrt(n)) noise level')
ylabel('Number of scans');
title('Noise reduction from averaging voxels')

%plot the level of noise for each scan averaged all voxels, so we can use for simulations
figure(156)
histogram(avgAllVoxelStdEachScan),
xlabel('Std of average of all voxel time series')
ylabel('Scans')
title('Std of time series averaged over all voxels, each scan')
keyboard

%% plot the snr curves
    %get the grad average of noise reduction as a function of voxel inclusion
    grandAvg = fittype('a/(x^b)', 'independent', 'x', 'dependent', 'y' );  
    gfit = fit(numAverages(numAverages<61)',voxelStdOverNumAverages(numAverages<61)',grandAvg)
    %plot it
    figure,scatter(numAverages,voxelStdOverNumAverages,1,'k','filled'), hold on
    plot(gfit),xlim([0 60]),ylim([0 3]);
    xlabel('Number of voxels')
    ylabel('Std of tseries')
    %plot snr curves
    figure, hold on
    expectedResponse = .169;
    for SNR = 1:10;
        fplot(@(x) ((gfit.a/(expectedResponse/SNR*sqrt(x)))^(1/gfit.b)))    
    end
    xlim([0 500]), ylim([0 300])
    xlabel('Number of scans')
    ylabel('Number of Voxels')
    legend('SNR 1','2','3','4','5')










%%%%%%%%%%
%% helpers%%
%%%%%%%%%%%

function [exponent, ratio, numAverages, voxelStdOverNumAverages] = checkAverageNoiseReduction(roi,avgSingleVoxelStd,subject,scan,numAverages,voxelStdOverNumAverages)

%average over different numbers of voxels and grab the standard deviations.
voxelAveragedStds = [];
%calculate the euclidian distance from the center
euclidDistFromCenter = sqrt((roi.scanCoords(1,:)-32).^2 + (roi.scanCoords(2,:)-32).^2);
[sortedDist distOrder] = sort(euclidDistFromCenter);
%average over voxels, adding them as a function of proximity
for numVoxels = 1:size(roi.tSeries,1);
    tSeries = roi.tSeries(distOrder(1:numVoxels),:);
    %average across all voxels
    tSeries = mean(tSeries,1)';
    %convert to % signal change
    meanTSeries = mean(tSeries);
    tSeries = 100*(tSeries-meanTSeries)/meanTSeries;
    voxelAveragedStds = [voxelAveragedStds std(tSeries)];
    %add to grand average
    numAverages = [numAverages numVoxels];
    voxelStdOverNumAverages = [voxelStdOverNumAverages std(tSeries)];
end
%set options
oneOverSqrtN = fittype('avgSingleVoxelStd/(x^b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0];
%fit a sqrt(n) improvement
oneOverSqrtNFit = fit((1:size(roi.tSeries,1))',voxelAveragedStds',oneOverSqrtN,opts);
% get the exponent of the fit - you would expect .5 if independent.
exponent = oneOverSqrtNFit.b;
%get the ratio of observed noise level to 1/sqrt(n) predicted noise level
ratio = voxelAveragedStds(end)/(avgSingleVoxelStd/sqrt(length(voxelAveragedStds)));

%plot a nice example
if subject == 5 & scan == 9;
    figure(155), subplot(1,2,1), hold on
    %plot std as a function of voxels averaged
    scatter(1:length(voxelAveragedStds),voxelAveragedStds,'MarkerFaceColor','k','MarkerEdgeColor','w');
    %plot the a/x^b fit
    plot(oneOverSqrtNFit);
    xlabel('Number of voxels averaged over')
    ylabel('Std of averaged time series')
    title('Example scan: std of numVoxels averaged')
    ylim([0 avgSingleVoxelStd]), xlim([0 length(voxelAveragedStds)+1])
    %plot actual 1/sqrt(x) improvement
    fplot(@(x) avgSingleVoxelStd./sqrt(x),'k')
    legend('Averages of different numVoxels','a/x^b fit','a/sqrt(x)')
end



