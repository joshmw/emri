function frequencyApproachSNR
cd '/Users/joshwilson/Library/CloudStorage/OneDrive-Stanford/emri/oldEmriSimulationData/simulData/frequencyAttempt'

%initialize empty average and SNRs
numSeeds = 150;
allScanTs = [];

%load all the time series
for seed = 1:numSeeds
    %load file and average time series over voxels
    fileName = sprintf('seed%4.4gfrequencyAttempt',seed); fileName = replace(fileName, '.', '_'); fileName = replace(fileName, ' ', '');
    load(fileName);
    brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);
    avgVoxelTs = squeeze(mean(ts(brain,brain,:),[1 2]));
    allScanTs(seed,:) = squeeze(avgVoxelTs);
end

%do some averaging
SNRs = [];
numSamples = 100;
for numScans = 1:numSeeds
    %loop through samples of different scans get an average of the average
    SNR = []
    for sample = 1:numSamples
        avgScanTs = mean(allScanTs(randperm(numSeeds,numScans),:),1);
        %calculate SNR from frequency components
        ft = abs(fft(avgScanTs/numScans));
        ft = ft(2:end);
        SNR = [SNR abs(ft(params.hz))/(abs(ft(params.hz-1))+abs(ft(params.hz+1))+abs(ft(params.hz-2))+abs(ft(params.hz+2)))*4];
    end
    %take the average average for that number of scans
    SNRs = [SNRs mean(SNR)];
    %graph if on the right iteration
    if mod(numScans,10) == 1;
        %graph time series
        figure(160), subplot(3,5,ceil(numScans/10))
        plot(squeeze(avgScanTs/numScans))
        ylim([99.5 100.5])
        title(sprintf('%i',numScans))
        set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[]), sgtitle('Time series over numAverages'),
        %graph frequency spectrum
        figure(161),  subplot(3,5,ceil(numScans/10))
        plot(abs(ft(1:round(end/2)))), hold on
        scatter(params.hz,ft(params.hz),'k','filled')
        set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
        title(sprintf('%i',numScans)), sgtitle('Frequency spectrums')
    end
end


%plot SNR change over number of averages
figure,scatter(1:150, SNRs)
xlabel('Number of scans averaged)')
ylabel('SNR (magnitude/adjacent 4 magnitudes)')
title('SNR over number of scans');


keyboard