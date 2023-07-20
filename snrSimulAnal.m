function snrSimulAnal

%go to data directory
dataDir = '/Users/joshwilson/Library/CloudStorage/OneDrive-Stanford/emri/oldEmriSimulationData/simulData/pinkNoiseSweep';
cd(dataDir)


%% get the base signal time series without noise
load('baseSignal.mat');
%get the brain and average the time series over it
brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);
ts = squeeze(ts);
noiselessTimeSeries = squeeze(mean(ts(brain,brain,:),[1 2]));

%% compare at different noise levels
figure(115), sub = 1;
for seed = 100
    permutedCorValsHigh = [];
    permutedCorValsLow = [];
    noiseStds = [];
    for pinkNoiseStd = [.025:.025:1];
        %load the file and squeeze+average the timeseries over the brain
        fileName = sprintf('seed%4.4gpinkNoise%4.4g',seed,pinkNoiseStd); fileName = replace(fileName, '.', '_'); fileName = replace(fileName, ' ', '');
        load(fileName);
        noisyTimeSeries = squeeze(ts);
        noisyTimeSeries = squeeze(mean(noisyTimeSeries(brain,brain,:),[1 2]));
        %plot it
        figure(115), sgtitle('Time series at different noise levels')
        subplot(5,8,sub),
        plot(noisyTimeSeries), ylim([99 102])
        title(sprintf('Std: %3.3g',pinkNoiseStd)), set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
        %plot the cross cor
        [cor lags] = xcorr(noisyTimeSeries-mean(noisyTimeSeries), noiselessTimeSeries-mean(noiselessTimeSeries),'Normalized');
        figure(116), sgtitle('Cross correlation with true signal at different noise levels')
        subplot(5,8,sub)
        plot(lags,cor), hold on, ylim([-.3 1])
        scatter(0,cor(200),'k','filled')
        title(sprintf('Std: %3.3g',pinkNoiseStd)), set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
        %plot the value? this is wrong
        figure(117), hold on,
        noiseStd = std(noisyTimeSeries(50:end)); noiseStds = [noiseStds noiseStd];
        scatter(noiseStd,corr2(noiselessTimeSeries,noisyTimeSeries),'k','filled');
        %permutation
        permutedCorVals = [];
        for perm = 1:500;
            shuffledNoisyTimeSeries = noisyTimeSeries(randperm(length(noisyTimeSeries)));
            permutedCorVals = [permutedCorVals corr2(noiselessTimeSeries,shuffledNoisyTimeSeries)];
        end
        permutedCorVals = sort(permutedCorVals);
        permutedCorValsHigh = [permutedCorValsHigh permutedCorVals(475)];
        permutedCorValsLow = [permutedCorValsLow permutedCorVals(25)];
        sub = sub+1;
    end
end

%get chance crosscor
[noiseStds, sortOrder] = sort(noiseStds);
plot(noiseStds,permutedCorValsHigh(sortOrder))
plot(noiseStds,permutedCorValsLow(sortOrder))
xlabel('Std of average time series'), ylabel('Correlation with; noiseless signal')
labelGraphs

%plot the injected vs reconstructed std
figure,scatter([1:40]/40,noiseStds)
xlabel('Std of injected pink noise in simulation')
ylabel('Std of averaged voxel time series (noise part)')


%% plot different number of averages for noisiest condition
sub = 1;
avgts = zeros(1,200);
noiseStdsAveraged = [];
for seed = 100:124
    for pinkNoiseStd = 1;
        %load the file and squeeze+average the timeseries over the brain
        fileName = sprintf('seed%4.4gpinkNoise%4.4g',seed,pinkNoiseStd); fileName = replace(fileName, '.', '_'); fileName = replace(fileName, ' ', '');
        load(fileName);
        noisyTimeSeries = squeeze(ts);
        noisyTimeSeries = squeeze(mean(noisyTimeSeries(brain,brain,:),[1 2]));
        %average the different runs and plot
        avgts = squeeze(avgts) + squeeze(noisyTimeSeries')./25;
        noiseStdsAveraged = [noiseStdsAveraged std(avgts(50:200).*(25/(seed-99)))];
        figure(118)
        subplot(5,5,sub)
        plot(avgts.*(25/(seed-99))), ylim([99 102]);
        sub = sub+1;
        %plot the correlation
        figure(119), hold on,
        scatter(seed-99,std(avgts(50:200).*(25/(seed-99)))),
        xlabel('Number of scans averaged'),
        ylabel('SNR')
    end
end
xlabel('Number of averaged time series'), ylabel('Correlation with noiseless signal')

keyboard








function labelGraphs
figure(115),subplot(5,8,17),
set(gca,'YTickLabel',[-1 0 1 2]), ylabel('% signal change')
subplot(5,8,36), 
set(gca,'XTickLabel',[0 500 1000]), xlabel('Time(ms)')

figure(116),subplot(5,8,17),
set(gca,'YTickLabel',[-.3 1]), ylabel('% correlation')
subplot(5,8,36), 
set(gca,'XTickLabel',[-200 200]), xlabel('Lag (ms)')







function getStdOfRealData

cd ~/data/emri/paperdata
SIDs = ["s0615" "s0616" "s0617" "s0618" "s0619"];

allVoxelStds = [];
avgStds = [];
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

        tSeries = roi.tSeries
        %average across all voxels
        tSeries = mean(tSeries,1)';
        %convert to % signal change
        meanTSeries = mean(tSeries);
        tSeries = 100*(tSeries-meanTSeries)/meanTSeries;
    
        %save the std of the averaged voxels in this scan
        avgStds = [avgStds std(tSeries)];

    end
    %go back and do it over
    pause(1)
    deleteView(v)
    mrQuit
    cd ..
   
end


function checkAverageNoiseReduction(roi)

voxelAveragedStds = []
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
oneOverSqrtN = fittype('a/sqrt(x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0];

%fit a sqrt(n) improvement
oneOverSqrtNFit = fit(1:size(roi.tSeries,1)',tSeries',oneOverSqrtN,opts);






