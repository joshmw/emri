function frequencyApproachSNR
% Using simultions, fits a numScans -> SNR function and inverts that along with the
% measured level of voxel noise to get isocontour SNR lines based on number of voxels and
% number of scans at a given frequency and response strength.
%
% USAGE:
%   %frequencyApproachSNR

%go to data folder
cd '/Users/joshwilson/Library/CloudStorage/OneDrive-Stanford/emri/oldEmriSimulationData/simulData/frequencyAttempt'
    
%initialize empty average and SNRs; set frequency you want to look at and number of simulations you have
frequency = 5;
numSeeds = 100;
allScanTs = [];

%load all the time series
for seed = 1:numSeeds
    %load file and average time series over voxels
    fileName = sprintf('Freq%iSeed%4.4gfrequencyAttempt',frequency,seed); fileName = replace(fileName, '.', '_'); fileName = replace(fileName, ' ', '');
    load(fileName);
    brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);
    avgVoxelTs = squeeze(mean(ts(brain,brain,:),[1 2]));
    allScanTs(seed,:) = squeeze(avgVoxelTs);
end

%do the averaging over different numbers of scans
SNRs = [];
numSamples = round(numSeeds*.5);
for numScans = 1:numSeeds
    %loop through samples of different scans get an average of the average
    SNR = [];
    for sample = 1:numSamples
        avgScanTs = mean(allScanTs(randperm(numSeeds,numScans),:),1);
        %calculate SNR from frequency components
        ft = abs(fft(avgScanTs/numScans));
        ft = ft(2:end);
        SNR = [SNR (abs(ft(params.hz)) - ((abs(ft(params.hz-1))+abs(ft(params.hz+1))+abs(ft(params.hz-2))+abs(ft(params.hz+2))))/4) / ...
            (abs(ft(params.hz-1))+abs(ft(params.hz+1))+abs(ft(params.hz-2))+abs(ft(params.hz+2)))*4];
        %SNR = [SNR (abs(ft(params.hz)) - ((abs(ft(params.hz-1))+abs(ft(params.hz+1))))/2) / ...
        %    (abs(ft(params.hz-1))+abs(ft(params.hz+1)))*2];
    end
    %take the average average for that number of scans
    SNRs = [SNRs mean(SNR)];
    %graph if on the right iteration
    if mod(numScans,10) == 1;
        %graph time series
        figure(160), subplot(3,5,ceil(numScans/10))
        plot(squeeze(avgScanTs))
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
figure(162), hold on
scatter(1:numSeeds, SNRs)
xlabel('Number of scans averaged)')
ylabel('SNR (magnitude/adjacent 4 magnitudes)')
title('SNR over number of scans');
%fit and plot a 1/sqrt(n)
oneOverSqrtN = fittype('a*sqrt(x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0];
oneOverSqrtNFit = fit((1:numSeeds)',SNRs',oneOverSqrtN,opts)
plot(oneOverSqrtNFit);

%plot SNR curve from the SNR fits
figure, hold on
SNR = 1; responseStrength = .169;
for SNR = 1:10
    fplot(@(x) ((1.706/(oneOverSqrtNFit.a*(responseStrength/.169)/SNR*sqrt(x)))^(1/.2263)))
    ylim([0 200]); xlim([0 300])
end


keyboard