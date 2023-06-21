function emriFrequencySpectrumAnal

%get the scan parameters
v = getMLRView;
ts = squeeze(loadTSeries(v));    

%get the voxels that have a high autocorrelation - last overlay of emriAnal
highAC = v.analyses{1}.overlays(end).data{v.curScan} > .5;

%get number of TRs, TR length, and x axis values of cycles/second
lengthInTRs = length(ts(1,1,:));
TRlength = viewGet(v,'framePeriod',v.curScan)
secondsInScan = lengthInTRs*TRlength;
freqComponents = 1:lengthInTRs-1;
cyclesPerSecond = freqComponents/secondsInScan;


%average the frequency components of all the voxels that meet the autocorrelation cutoff
allFtSeries = zeros(1,round(lengthInTRs/2)-1);
for row = 1:64 
    for col = 1:64
        if highAC(row,col) == 1;
            ftSeries = squeeze(abs(fft(ts(row,col,:))));
            ftSeries = ftSeries(2:round(lengthInTRs/2));
            ftSeries = ftSeries/sum(ftSeries);
            allFtSeries = allFtSeries + ftSeries;
        end
    end
end

%pick a line type based on the name of the scan you are looking at - gated vs ungated vs flipped encoding direction are different
if contains(v.figure.Name,'cardiacGated')
    lineType = '-';
elseif contains(v.figure.Name,'cardiacUngated')
    lineType = '--';
elseif contains(v.figure.Name,'flipped');
    lineType = '-.'
else
    sprintf('Not sure which scan this is. Description should include "cardiacGated", "cardiacUngated", or "flipped".')
    keyboard
end

%plot the power vs the cycle length
plot(cyclesPerSecond(1:length(allFtSeries)),allFtSeries/sum(sum(highAC)),lineType)
xlabel('Cycles/second')
ylabel('Normalized power');
