function [frequencies, amplitudes] = emriFrequencySpectrumAnal(v,isSimul,color,frequencies,amplitudes,highACVal)
% to run as a standalone on the current scan in your view:
% v = getMLRView;
% emriFrequencySpectrumAnal([],1,'k',[],[],.5) 


%get the view if you didnt pass one in
if isempty(v);
    v = getMLRView;
end

%load the time series
ts = squeeze(loadTSeries(v));    

%get the voxels that have a high autocorrelation - last overlay of emriAnal
highAC = v.analyses{1}.overlays(end).data{v.curScan} > highACVal;
%get number of TRs, TR length, and x axis values of cycles/second
lengthInTRs = length(ts(1,1,:));
TRlength = viewGet(v,'framePeriod',v.curScan);
if isSimul;
    TRlength = .005;
end
secondsInScan = lengthInTRs*TRlength;
freqComponents = 1:lengthInTRs-1;
cyclesPerSecond = freqComponents/secondsInScan;

%set the number of frequency components to look at - this accounts for different length scans
numComponents = 50;
if numComponents > lengthInTRs/2;
    sprintf('The number of frequency components you are plotting is high that what you have in the actual data.')
    keyboard
end

%average the frequency components of all the voxels that meet the autocorrelation cutoff
allFtSeries = zeros(1,numComponents);
for row = 1:length(highAC)
    for col = 1:length(highAC)
        if highAC(row,col) == 1
            ftSeries = squeeze(abs(fft(ts(row,col,:))))';
            ftSeries = ftSeries(2:numComponents+1);
            ftSeries = ftSeries/sum(ftSeries);
            allFtSeries = allFtSeries + ftSeries;
        end
    end
end

%pick a line type based on the name of the scan you are looking at - gated vs ungated vs flipped encoding direction are different
lineType = '-';

%plot the power vs the cycle length
amps = allFtSeries/sum(sum(highAC));
if ~isSimul, figure(150), hold on,
xlabel('Frequency (cycles/second)'), ylabel('Normalized Amplitude'),
plot(cyclesPerSecond(1:numComponents),amps,'lineStyle',lineType,'Color',color);
figure(151), hold on, xlabel('Frequency (cycles/second)'), ylabel('Normalized Power');
end

plot(cyclesPerSecond(1:numComponents),amps.^2,'lineStyle',lineType,'Color',color);

%ylabel('Normalized power');
frequencies = [frequencies cyclesPerSecond(1:numComponents)];
amplitudes = [amplitudes amps];



