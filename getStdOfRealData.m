function getStdOfRealData

cd ~/data/emri/paperdata
SIDs = ["s0615" "s0616" "s0617" "s0618" "s0619"];

allVoxelStds = [];
avgAllVoxelStdEachScan = [];
voxelAverageImprovementFitsEachScan = [];

%get all of the data
for subject = 1:length(SIDs);

    cd(SIDs(subject));
    %create the view and go to unaveraged ungated scana
    v = newView;
    v = viewSet(v,'curGroup','rawUngated')
    v = loadROI(v,'insideBrain') ;
    numScans = viewGet(v,'numScans');
    %loop through each scan
    for scan = 1:numScans
        %set the current scan
        v = viewSet(v,'curScan',scan)    
        %load the t series
        roi = loadROITSeries(v,'insideBrain',viewGet(v,'curScan'),viewGet(v,'curGroup'));

        %invidiaul t series std
        alltSeriesMeans = mean(roi.tSeries,2);
        alltSeries = 100*(roi.tSeries-alltSeriesMeans)./repmat(alltSeriesMeans,[1 size(roi.tSeries,2)]);
        allVoxelStds = [allVoxelStds std(alltSeries,0,2)'];

        %average across all voxels
        tSeries = roi.tSeries;
        tSeries = mean(tSeries,1)';
        %convert to % signal change
        meanTSeries = mean(tSeries);
        tSeries = 100*(tSeries-meanTSeries)/meanTSeries;
        %save the std of the averaged voxels in this scan
        avgAllVoxelStdEachScan = [avgAllVoxelStdEachScan std(tSeries)];

        %check the improvement averaging over each scan
        exponent = checkAverageNoiseReduction(roi, subject, scan);
        voxelAverageImprovementFitsEachScan = [voxelAverageImprovementFitsEachScan exponent];

    end
    %go back and do it over
    pause(1)
    deleteView(v)
    mrQuit
    cd ..
   
end

%plot the improvement with averaging more voxels
figure(155), subplot(1,2,2)
hist(voxelAverageImprovementFitsEachScan); xlim([0 1])
xlabel('Improvement parameter: b in (a/(x^b))')
ylabel('Number of scans');
title('Noise reduction from averaging voxels')

%plot the level of noise for each scan averaged all voxels, so we can use for simulations
figure(156)
hist(avgAllVoxelStdEachScan),
xlabel('Std of average of all voxel time series')
ylabel('Scans')
title('Std of time series averaged over all voxels, each scan')
keyboard



function exponent = checkAverageNoiseReduction(roi,subject,scan)

%average over different numbers of voxels and grab the standard deviations.
voxelAveragedStds = [];
for numVoxels = 1:size(roi.tSeries,1);
    tSeries = roi.tSeries(randperm(size(roi.tSeries,1),numVoxels),:);
    %average across all voxels
    tSeries = mean(tSeries,1)';
    %convert to % signal change
    meanTSeries = mean(tSeries);
    tSeries = 100*(tSeries-meanTSeries)/meanTSeries;
    voxelAveragedStds = [voxelAveragedStds std(tSeries)];
end

%set options
oneOverSqrtN = fittype('a/(x^b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0];
%fit a sqrt(n) improvement
oneOverSqrtNFit = fit((1:size(roi.tSeries,1))',voxelAveragedStds',oneOverSqrtN,opts);
% get the exponent of the fit - you would expect .5 if independent.
exponent = oneOverSqrtNFit.b;


%plot a nice example
if subject == 5 & scan == 9;
    figure(155), subplot(1,2,1), hold on
    scatter(1:length(voxelAveragedStds),voxelAveragedStds,1,'k');
    plot(oneOverSqrtNFit);
    xlabel('Number of voxels averaged over')
    ylabel('Std of averaged time series')
    title('Example scan: std of numVoxels averaged')
    ylim([0 1]), xlim([0 250])
    legend('Averages of different numVoxels','a/x^b fit')
    figure
end



