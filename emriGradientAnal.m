function [brainMean, spreadMean, outsideMean, brainStd, spreadStd, outsideStd, brainNum, spreadNum, outsideNum, brainAllAutocor, spreadAllAutocor, outsideAllAutocor] = emriGradientAnal(v)

sprintf('Only do this one on real data!! Use emriSimulAnal for simulation data. This will pick up on ROIs drawn on real data.')

%set a few parameters that you probably don't want to change, but can
graphStuff = 1;
frequencyToView = 1;

%get information about the view
scan = v.curScan
if scan == 3;
    rois = [4 5 6]
else
    rois = [1 2 3]
end
params.xdim = 64;

%check you have the right ROIs
if ~(contains(v.ROIs(rois(1)).name,'spread') & contains(v.ROIs(rois(2)).name,'insideBrain') & contains(v.ROIs(rois(3)).name,'outsideBrain'))
    sprintf('Unexpected ROIs. ROIs should be (in order): 1: "spread": A rectangle of spread around the brain in the phase-encoding direction; 2: "brain" - the brain in the slice; ', ...
        '3: "outsideBrain" - spread minus brain. So the voxels that have artifacts but are not in the brain. You do not need to define the voxels that are outside brain and not in spread direction')
    keyboard
end

%get the scan coords of the voxels in the roi
insideBrain = loadROITSeries(v,v.ROIs(rois(2)).name);
insideBrainCoords = insideBrain.scanCoords;
outsideBrain = loadROITSeries(v,v.ROIs(rois(3)).name);
outsideBrainCoords = outsideBrain.scanCoords;

%make masks for overlays
brainMask = zeros(64,64); spreadMask = zeros(64,64); outsideMask = zeros(64,64);

    %get the brain voxels
    for voxel = 1:length(insideBrainCoords)
        brainMask(insideBrainCoords(1,voxel,1,1), insideBrainCoords(2,voxel,1,1)) = 1;
    end
    
    %get the voxels along the spread dimension
    for voxel = 1:length(outsideBrainCoords);
        spreadMask(outsideBrainCoords(1,voxel), outsideBrainCoords(2,voxel)) = 1;
    end

    %get the voxels outside of brain and spread
    outsideMask = outsideMask - brainMask - spreadMask + 1;


%% vector average and fit a cosine to the phase gradient %%
graphGradient = 1;
if graphGradient
    %figure out the brain area to vector average over
    brain = min(insideBrainCoords(2,:)):max(insideBrainCoords(2,:));
    
    %average and plot
    figure, meanAmps = []; meanPhases = []; ampMeanPhases = [];
    %loop through relevant voxels by row - only works for x encoding direction
    for row = 1:params.xdim
        vectorsX = []; vectorsY = []; amplitudes = []; ampVectorsX = []; ampVectorsY = [];
        for col = brain
            %grab the amplitude and phase
            amplitude = v.analyses{1}.overlays(frequencyToView*3-1).data{scan}(row,col);
            phase = v.analyses{1}.overlays(frequencyToView*3).data{scan}(row,col);
            
            %grab the individual vectors, with and without amplitude
            vectorsX = [vectorsX cos(phase)];
            vectorsY = [vectorsY sin(phase)];
            ampVectorsX = [ampVectorsX amplitude*cos(phase)];
            ampVectorsY = [ampVectorsY amplitude*sin(phase)];    
            amplitudes = [amplitudes amplitude];
    
        end
            %average the vectors, with and without amplitude
            meanX = mean(vectorsX);
            meanY = mean(vectorsY);
            ampMeanX = mean(ampVectorsX);
            ampMeanY = mean(ampVectorsY);
            meanAmplitude = sqrt(meanX^2+meanY^2);
            meanPhase = cart2pol(meanX,meanY);
            ampMeanPhase = cart2pol(ampMeanX,ampMeanY);
    
            %plot the average vector in polar space and get the sequential averaged phase
            subplot(2,12,4:7)
            polarplot([0 meanPhase], [0 meanAmplitude]), title('Vector averages across rows'), hold on, rlim([0 1])
            meanAmps = [meanAmps meanAmplitude];
            meanPhases = [meanPhases meanPhase];
            ampMeanPhases = [ampMeanPhases ampMeanPhase];
    end
    
        %plot the phase map
        subplot(2,12,1:2),
        imagesc(v.analyses{1}.overlays(frequencyToView*3).data{scan})
        title('Phase map'); xlabel('x'), ylabel('y'),colorbar, colormap(jet), caxis([0 2*pi])
    
        %plot the vector averaged row phase map
        subplot(2,12,3)
        subZeroPhases = meanPhases(:)<0;
        imagesc(meanPhases(:)+subZeroPhases*2*pi)
        title('VA phase map'), xlabel('x'), ylabel('y'), colorbar, colormap(jet), caxis([0 2*pi])
    
        %plot a histogram of the averaged amplitudes 
        subplot(2,12,[8 9])
        histogram(meanAmps,'binEdges',[0.001:.1:1.001]), xlim([0 1])
        xlabel('Amplitude'), ylabel('count'), title('Amplitude')
        
        %plot a histogram of the averaged phases
        subplot(2,12,[11 12]);
        histogram(meanPhases,'binEdges',[-pi:pi/10:pi]), xlim([-pi-.1 pi+.1])
        xlabel('Phase (rad)'), title('Phase'); 
        
    
        %quantify gradiant with amplitude of best fitting cosine function
            %make y
            y = ampMeanPhases;
            y = y(:);
            %make x go between 0 and 2pi, non-inclusive
            x = linspace(0,2*pi,length(y)+1);
            x = x(1:end-1);
            x = x(:);
            %find the best frequency to describe the gradient by maximizing r2
            testFrequencies = .05:.05:5; r2s = [];
            for harmonicFreq = testFrequencies
                A = ones(length(y),3);
                A(:,1) = cos(harmonicFreq*x);
                A(:,2) = sin(harmonicFreq*x);
                %compute least squares estimate for harmonic frequency
                estimate = inv(A'*A)*A'*cos(y);
                %get r2
                residualVariance = var(cos(y)-A*estimate);
                originalVariance = var(cos(y));
                r2 = 1 - residualVariance/originalVariance;
                r2s = [r2s r2];
                %get theta and amplitude
                [theta r] = cart2pol(estimate(1),estimate(2));
                theta = rad2deg(theta);
                cosAmplitude = estimate(2);
                %plot r2
                subplot(2,12,[17 18]), hold on,
                scatter(harmonicFreq,r2), ylim([0 1]); xlabel('Model Cosine frequency'),ylabel('Model r2'),title('r2 over cos frequency sweep')
                %plot gradient cosine amplitude
                subplot(2,12,[ 20 21]), hold on,
                scatter(harmonicFreq,cosAmplitude), ylim([-1 1]), xlabel('Model cosine frequency'),ylabel('Model Amplitude'),title('Amp over cos frequency sweep')
                %plot gradient phase
                subplot(2,12,[ 23 24]), hold on,
                scatter(harmonicFreq,theta), xlabel('Model cosine frequency'),ylabel('Model cosine phase (def)'),title('Phase over cos frequency sweep')
            end
            
            %grab the best frequency and get the fit
            [maxr2, bestFrequencyIndex] = max(r2s),
            bestFrequency = testFrequencies(bestFrequencyIndex),
            A = ones(length(y),3);
            A(:,1) = cos(bestFrequency*x);
            A(:,2) = sin(bestFrequency*x);
            estimate = inv(A'*A)*A'*cos(y);
            [theta r] = cart2pol(estimate(1),estimate(2));
            bestTheta = rad2deg(theta);
            bestCosAmplitude = estimate(2);
            %plot the gradient and best cosine fit
            subplot(2,12,[13:15]), hold on
            plot(cos(y),'Color','b')
            plot(1:length(y),A*estimate,'--','Color','b');
            ylim([-1 1])
            xlabel('Sequential voxel'),
            ylabel('Cosine of vector-averaged phase');
            title('Phase gradient along spread direction')
            legend({'cosine','model'});
end

%% do the overlay analyses - make histograms %%

   %grab the overlays
    amplitudeOverlay = v.analyses{1}.overlays(frequencyToView*3-1).data{scan};
    phaseOverlay = v.analyses{1}.overlays(frequencyToView*3).data{scan};
    coherenceOverlay = v.analyses{1}.overlays(frequencyToView*3-2).data{scan};
    autocorOverlay = v.analyses{1}.overlays(end).data{scan};

    %plot the histograms for different sections
    %inside brain
    figure,
    subplot(2,2,1),histogram(amplitudeOverlay(brainMask>0),'FaceColor','g','BinEdges',[0:.2:5],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 5]), hold on
    xlabel('Amplitude'),title('Amplitude overlay distribution'), ylabel('Density')
    subplot(2,2,2),histogram(coherenceOverlay(brainMask>0),'FaceColor','g','BinEdges',[0:.05:1],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 1]), hold on
    xlabel('Coherence'),title('Coherence overlay distribution'), ylabel('Density')
    subplot(2,2,3),histogram(phaseOverlay(brainMask>0),'FaceColor','g','BinEdges',[0:pi/10:2*pi],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 2*pi]), hold on
    xlabel('Phase'),title('Phase overlay distribution'), ylabel('Density')
    subplot(2,2,4),histogram(autocorOverlay(brainMask>0),'FaceColor','g','BinEdges',[0:.05:1],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 1]), hold on
    xlabel('Autocorrelation'),title('Autocorrelation overlay distribution'), ylabel('Density')
    %along spread     
    subplot(2,2,1),histogram(amplitudeOverlay(spreadMask>0),'FaceColor','c','BinEdges',[0:.2:5],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 5])
    subplot(2,2,2),histogram(coherenceOverlay(spreadMask>0),'FaceColor','c','BinEdges',[0:.05:1],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 1])
    subplot(2,2,3),histogram(phaseOverlay(spreadMask>0),'FaceColor','c','Binedges',[0:pi/10:2*pi],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 2*pi])
    subplot(2,2,4),histogram(autocorOverlay(spreadMask>0),'FaceColor','c','BinEdges',[0:.05:1],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 1])
    %outside
    subplot(2,2,1),histogram(amplitudeOverlay(outsideMask>0),'FaceColor','m','BinEdges',[0:.2:5],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 5])
    subplot(2,2,2),histogram(coherenceOverlay(outsideMask>0),'FaceColor','m','BinEdges',[0:.05:1],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 1])
    subplot(2,2,3),histogram(phaseOverlay(outsideMask>0),'FaceColor','m','BinEdges',[0:pi/10:2*pi],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 2*pi])
    subplot(2,2,4),histogram(autocorOverlay(outsideMask>0),'FaceColor','m','BinEdges',[0:.05:1],'FaceAlpha',.3,'EdgeColor','none','Normalization','probability'); xlim([0 1])
   
%get the summary values for emriGradientAnalPlot
brainMean = mean(autocorOverlay(brainMask>0)); brainStd = std(autocorOverlay(brainMask>0)), brainNum = sum(sum(brainMask>0));
spreadMean = mean(autocorOverlay(spreadMask>0)); spreadStd = std(autocorOverlay(spreadMask>0)), spreadNum = sum(sum(spreadMask>0));
outsideMean = mean(autocorOverlay(outsideMask>0)); outsideStd = std(autocorOverlay(outsideMask>0)), outsideNum = sum(sum(outsideMask>0));

%pull out all autocorData
brainAllAutocor = autocorOverlay(brainMask>0);
spreadAllAutocor = autocorOverlay(spreadMask>0);
outsideAllAutocor = autocorOverlay(outsideMask>0);

%get the p values for the spread to say it's in one direction
circStdsAlongSpread = []
for line = brain;
    circStdsAlongSpread = [circStdsAlongSpread circ_std(phaseOverlay(:,line))];
end

circStdsAgainstSpread = []

for line = 1:64;
    circStdsAgainstSpread = [circStdsAgainstSpread circ_std(phaseOverlay(line,brain)')];
end

figure, subplot(1,2,1), hist(circStdsAlongSpread),subplot(1,2,2), hist(circStdsAgainstSpread)

[h p] = ttest2(circStdsAgainstSpread,circStdsAlongSpread,'tail','left') 
title(sprintf('p = %4.4i',p))






%% helpers %%
function plotSummary

figure(110), hold on
violin(coherenceOverlay(brainMask>0),'facecolor','g'), ylim([0 1])
violin(coherenceOverlay(spreadMask>0),'facecolor','c'), ylim([0 1])
violin(coherenceOverlay(outsideMask>0),'facecolor','m'), ylim([0 1]), painters

figure(111), hold on
violin(autocorOverlay(brainMask>0),'facecolor','g'), ylim([0 1])
violin(autocorOverlay(spreadMask>0),'facecolor','c'), ylim([0 1])
violin(autocorOverlay(outsideMask>0),'facecolor','m'), ylim([0 1]), painters

figure(112), hold on
violin(phaseOverlay(brainMask>0),'facecolor','g'), ylim([0 2*pi])
violin(phaseOverlay(spreadMask>0),'facecolor','c'), ylim([0 2*pi])
violin(phaseOverlay(outsideMask>0),'facecolor','m'), ylim([0 2*pi]), painters

sid = 's0619'; gating = 'ungated';
savepdf(figure(110),strcat(sid,gating,'CoherenceSummary'))
savepdf(figure(111),strcat(sid,gating,'AutocorSummary'))
savepdf(figure(112),strcat(sid,gating,'PhaseSummary'))


%scatter plots

cohBrain = mean(coherenceOverlay(brainMask>0))
cohSpread = mean(coherenceOverlay(spreadMask>0))
cohOutside = mean(coherenceOverlay(outsideMask>0))

autocorBrain = mean(autocorOverlay(brainMask>0))
autocorSpread = mean(autocorOverlay(spreadMask>0))
autocorOutside = mean(autocorOverlay(outsideMask>0))

subplot(1,2,1), hold on
title('Flipped')
scatter(1,autocorBrain,'MarkerFaceColor','g','MarkerEdgeColor','w')
scatter(2,autocorSpread,'MarkerFaceColor','c','MarkerEdgeColor','w')
scatter(3,autocorOutside,'MarkerFaceColor','m','MarkerEdgeColor','w')
xlim([0 4]), ylim([0 1])



subplot(1,2,1), hold on
title('Gated')
scatter(1,autocorBrain,'MarkerFaceColor','g','MarkerEdgeColor','w')
scatter(2,autocorSpread,'MarkerFaceColor','c','MarkerEdgeColor','w')
scatter(3,autocorOutside,'MarkerFaceColor','m','MarkerEdgeColor','w')
xlim([0 4]), ylim([0 1])

subplot(1,2,2), hold on
title('Ungated')
scatter(1,autocorBrain,'MarkerFaceColor','g','MarkerEdgeColor','w')
scatter(2,autocorSpread,'MarkerFaceColor','c','MarkerEdgeColor','w')
scatter(3,autocorOutside,'MarkerFaceColor','m','MarkerEdgeColor','w')
xlim([0 4]), ylim([0 1])



