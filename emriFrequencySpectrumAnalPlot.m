function [frequencies, amplitudes] = emriFrequencySpectrumAnalPlot(SIDs,highACVal)
% cd ~/data/emri/paperdata
%  emriFrequencySpectrumAnalPlot(["s0615" "s0616" "s0617" "s0618" "s0619"],.7)

cd ~/data/emri/paperdata
frequencies = []; 
amplitudes = [];

for subject = 1:length(SIDs);
    cd(SIDs(subject));
    %create the view and go to ungated scan
    v = newView;
    v = viewSet(v,'curGroup','averagesResized')
    v = viewSet(v,'curScan',2)
    v = loadAnalysis(v,'emriAnal');

    %get the frequencies and amplitudes
    [frequencies, amplitudes] = emriFrequencySpectrumAnal(v,0,'k',frequencies,amplitudes,highACVal);
    
    %go back and do it over
    pause(1)
    deleteView(v)
    mrQuit
    cd ..
   
end

% sort it
[frequencies, sortorder] = sort(frequencies);
amplitudes = amplitudes(sortorder)

%fit a 1/f spectrum
oneOverF = fittype('a*(1/x) + b', 'independent', 'x', 'dependent', 'y' );
oneOverFSquared = fittype('(a*(1/x^2) + b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0];

[ampFit ampFitGof] = fit(frequencies',amplitudes',oneOverF,opts)
[powerFit powerFitGof] = fit(frequencies',amplitudes.^2',oneOverFSquared,opts)

%show the amp fits
figure(150),
plot(ampFit,'r'),ylim([0 .12])
copy, xlim([0 15])

%plot the power fits
figure(151),
plot(powerFit,'r'), ylim([0 .015])
copy, xlim([0 15])

keyboard


%% plot the simulations
cd '~/data/EMRI/paperData/simulations/mrsesh'

nScans = 101;

lightBLUE = [0.356 0.811 0.956];
darkBLUE = [0.0196 0.074 0.670];
diffBLUE = lightBLUE - darkBLUE;

figure, hold on

    v = newView;
    v = viewSet(v,'curGroup','freqSweep');
    v = loadAnalysis(v,'emriAnal');


for scan = 2:4:99,
    v = viewSet(v,'curScan',scan); 


    subplot(5,5,ceil(scan/4)), hold on,

    emriFrequencySpectrumAnal(v,1,[darkBLUE+scan/nScans*diffBLUE],[],[],highACVal);

    plot(powerFit,'r');
    
    xlim([0 15]),ylim([0 1]),
    xticklabels([]), yticklabels([]);

    %drawPublishAxis, legend('off')

    pause(1)
end
