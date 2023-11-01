function plotEmriTimeSeries
cd '/Users/joshwilson/data/EMRI/paperData'
SIDs = ["s0615" "s0616" "s0617" "s0618" "s0619"];
figure(1), hold on

%get the data and plot raw time series.
for subject = 1:length(SIDs);

    cd(SIDs(subject))
    v = newView;
    v = viewSet(v,'curGroup','averagesResized')
    v = viewSet(v,'curScan',1)

    %load tseries
    v123 = loadROITSeries(v,'v123') 
    
    %convert to % signal change
    numVoxels = height(v123.tSeries);
    nSamples = width(v123.tSeries);
    t = (1:nSamples)*3.3;
    v123.tSeriesPercent = v123.tSeries;
    for vox = 1:numVoxels  
        v123.tSeriesPercent(vox,1:nSamples) = 100*(v123.tSeries(vox,:)-mean(v123.tSeries(vox,:)))/mean(v123.tSeries(vox,:));
    end
    %plot individual voxels in the ROI
    figure, for vox = 1:min(100,height(v123.tSeriesPercent)),subplot(10,10,vox),plot(t,squeeze(v123.tSeriesPercent(vox,:))),end
    sgtitle('V1-V3 time series, individual voxels, averaged scans')
    subplot(10,10,51),ylabel('% signal change'); subplot(10,10,95),xlabel('Time (ms)')
    %plot the average
    avgTSeries = squeeze(mean(v123.tSeriesPercent));
    figure,plot(t,avgTSeries)   
    title('Averaged scans, averaged voxels V1-V3')
    xlabel('Time (ms/3)')
    ylabel('% signal change')
    figure(1), plot(t,avgTSeries) 

    %save the time series and parameters of the scan
    allVoxelTSeries{subject} = v123.tSeriesPercent;
    allV123{subject} = v123;

    %move on to other subjecrts
    pause(1)
    deleteView(v)
    mrQuit
    cd ..
end

keyboard

% plot the frequency spectrum of all the voxels for each subject
figure, hold on
for subject = 1:length(SIDs)
    %initialize empty frequency spectrum
    ft = zeros(1,width(allVoxelTSeries{subject}));
    %add all the frequency spectrums from all the voxels
    for voxel = 1:height(allVoxelTSeries{subject})
        ft = ft + abs(fft(allVoxelTSeries{subject}(voxel,:)))
    end
    ft = ft(2:51);
    ft = ft/sum(ft);
    scanLength = allV123{subject}.framePeriod * width(allVoxelTSeries{subject});
    plot([1:50]/scanLength,ft)
end
xlim([0 20])
xlabel('Cycles per second')
ylabel('Amplitude')
title('Averaged frequency spectrum of all voxels')
legend('Subject 1','2','3','4','5')






function plotSVD
[U S V] = svd(x)
S2 = S;   
for i  = 4:65, S2(i,i) = 0; end 
Y2 = U * S2 * V';    
























% %show the average of all of the scans+voxels for all subjects
% load("allScanTs.mat");
% figure
% for subject = 1:5,
%     subplot(1,5,subject)
%     plot(mean(squeeze(allSubjectAvgTs{subject}),1))
%     ylim([-1 1])
%     xlabel('Time (ms)'),
%     ylabel('% signal change')
%     title(sprintf('Subject %i', subject))
% end


