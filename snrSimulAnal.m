function snrSimulAnal

showVoxelAveraging = 0;
showMeanAveraging = 1;


%go to data directory
dataDir = '/Users/joshwilson/Library/CloudStorage/OneDrive-Stanford/emri/oldEmriSimulationData/simulData/pinkNoiseSweep';
cd(dataDir)


%% get the base signal time series without noise
load('baseSignal.mat');
%get the brain and average the time series over it
brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);
ts = squeeze(ts);
noiselessTimeSeries = squeeze(mean(ts(brain,brain,:),[1 2]));

if showVoxelAveraging
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
        %sort the permuted values and get upper/lower bounds
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

end



%% show the effects of averaging on data simulated to look like human data
%go to data
dataDir = '/Users/joshwilson/Library/CloudStorage/OneDrive-Stanford/emri/oldEmriSimulationData/simulData/snrSimul';
cd(dataDir)

colors{1} = [0 208 233]/255;
colors{2} = [126 200 42]/255;
pinkNoiseStd = [[1.07 1.4]];

for noiseLevel = 1:length(pinkNoiseStd)
    avgts = zeros(1,200);
    noiseStdsAveraged = [];
    for seed = 1:200
        %load the file and squeeze+average the timeseries over the brain
        fileName = sprintf('seed%4.4gpinkNoise%4.4g',seed,pinkNoiseStd(noiseLevel)); fileName = replace(fileName, '.', '_'); fileName = replace(fileName, ' ', '');
        load(fileName);
        noisyTimeSeries = squeeze(ts);
        noisyTimeSeries = squeeze(mean(noisyTimeSeries(brain,brain,:),[1 2]));
        %average the different runs and plot
        avgts = squeeze(avgts) + squeeze(noisyTimeSeries')./200;
        noiseStdsAveraged = [noiseStdsAveraged std(avgts(50:200).*(200/(seed)))];
        %plot if
        if noiseLevel == 1;
            if (seed == 12) | (seed == 47) | (seed == 106) | (seed == 188)
                figure
                plot(avgts.*(200/(seed))), ylim([99.8 100.5]);
            end
        end
    end
    figure(125); hold on
    scatter(1:200,noiseStdsAveraged,'filled','MarkerFaceColor',colors{noiseLevel},'MarkerEdgeColor','w'),   
    title('Influence of averaging over scans on noise magnitude')
    xlabel('Number of scans averaged')
    ylabel('Std of time series outside of signal')
    %fit
    oneOverSqrtN = fittype('a/(sqrt(x))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); opts.Display = 'Off'; opts.StartPoint = [0];
    oneOverSqrtNFit = fit([1:length(noiseStdsAveraged)]',noiseStdsAveraged',oneOverSqrtN,opts);
    plot(oneOverSqrtNFit)
    xlim([0 100]); ylim([0 .75])
end
plot([0 100], [.169/2 .169/2])

keyboard


%%%%%%%%%%%%%
%% helpers %%
%%%%%%%%%%%%%

function labelGraphs
figure(115),subplot(5,8,17),
set(gca,'YTickLabel',[-1 0 1 2]), ylabel('% signal change')
subplot(5,8,36), 
set(gca,'XTickLabel',[0 500 1000]), xlabel('Time(ms)')

figure(116),subplot(5,8,17),
set(gca,'YTickLabel',[-.3 1]), ylabel('% correlation')
subplot(5,8,36), 
set(gca,'XTickLabel',[-200 200]), xlabel('Lag (ms)')



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






