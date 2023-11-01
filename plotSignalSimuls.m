function plotSignalSimuls(seeds,offsets)
% plotSignalSimuls(100:120,1:25)

%go to data
cd '~/Library/CloudStorage/OneDrive-Stanford/emri/oldEmriSimulationData/simulData/signalTesting'

%initiate empty matrices of fwhm and max peaks
nFrames = 40; %hard coded
maxPeakValsAll = [];
fwhmValsAll = [];
peakLocationsAll = [];
avgTSeries = zeros(length(offsets),nFrames);

%sweep through the different randomizations
for seed = seeds

    %initialize empty arrays
    maxPeakVals = [];
    fwhmVals = [];

    %sweep through the signal offset magnitudes
    for offset = offsets
        
        %load the simulation
        fileName = sprintf('seed%4.4gOffset%4.4g',seed,offset);
        fileName = replace(fileName, '.', '_');
        fileName = replace(fileName, ' ', '');
        load(fileName);

        %brain voxels
        brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);
        
        %loop through voxels and get the average peak and FWHM
        maxPeakValsSingleOffset = [];
        fwhmValsSingleOffset = [];
        for voxelx = brain,
            for voxely = brain;
                %average every voxel in the brain
                tsVoxel = squeeze(ts(voxelx,voxely,1,:)) - params.brainBaseContrast;
                maxPeakValsSingleOffset = [maxPeakValsSingleOffset max(tsVoxel)];
                fwhmValsSingleOffset = [fwhmValsSingleOffset findFWHM(1:length(tsVoxel),tsVoxel)];
                [maxVal maxLoc] = max(tsVoxel);
                peakLocations{offset}(voxelx, voxely) = maxLoc;
                %get the average time series of every voxel/seed in the brain for each offset
                avgTSeries(offset,:) = avgTSeries(offset,:) + tsVoxel'./(length(brain)^2*length(seeds));
            end
        end
        %add the average peak/fwhm of all voxels for that offset to the data
        maxPeakVals = [maxPeakVals mean(maxPeakValsSingleOffset)];
        fwhmVals = [fwhmVals mean(fwhmValsSingleOffset)];
    end

    maxPeakValsAll(seed-99,:) = maxPeakVals;
    fwhmValsAll(seed-99,:) = fwhmVals;

    %plot the peak values across offsets
    figure(1), hold on,
    scatter(offsets,maxPeakVals,'filled','MarkerEdgeColor','w')
    xlabel('Signal randomization (ms)'), ylabel('Signal peak value')
    %plot the fwhm values across offsets
    figure(2), hold on
    scatter(offsets,fwhmVals*params.sampleTimeMs,'filled','MarkerEdgeColor','w')
    xlabel('Signal randomization (ms)'), ylabel('Signal FWHM max (ms)')

end

%plot summary lines
figure(1), painters,
plot(1:25,mean(maxPeakValsAll,1),'Color','k','MarkerEdgeColor','w')

figure(2), painters,
plot(1:25,mean(fwhmValsAll,1)*params.sampleTimeMs,'Color','k','MarkerEdgeColor','w')



%% plot the actual time series
figure, hold on

%reinterpolating back up to original signal x dimension (ms)
t = linspace(0, 1, params.numKsamples);
n = params.ntimePointsMs*2; % required length
t_ = linspace(0, 1, n);

% plot with interpolation + shift so they are all centered. need to interpolated to 2x the original signal resolution so you can resolve
%the ambiguity of the signal being centered between time points (i.e. 101.5/200ms -> 203/400ms)
for offset = 1:25,
    color = offset/(35);
    plot(circshift(interp1(t, avgTSeries(offset,:), t_, 'spline').',-offset),'Color',[color color color]),
end

keyboard
