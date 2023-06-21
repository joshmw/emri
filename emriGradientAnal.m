function emriGradientAnal

%set a few parameters that you probably don't want to change, but can
graphStuff = 1;
frequencyToView = 1;

%get information about the view
v = getMLRView;
scan = v.curScan
params.xdim = length(v.baseVolumes(v.curBase).data);

%figure out the brain area to vector average over
brain = min(v.ROIs(1).coords(2,:)) : max(v.ROIs(1).coords(2,:));

%make masks for overlays
brainMask = zeros(64,64); spreadMask = zeros(64,64); outsideMask = zeros(64,64);
    
    %get the brain voxels
    for voxel = 1:length(v.ROIs(2).coords);
        brainMask(v.ROIs(2).coords(1,voxel), v.ROIs(2).coords(2,voxel)) = 1;
    end
    
    %get the voxels along the spread dimension
    for voxel = 1:length(v.ROIs(3).coords);
        spreadMask(v.ROIs(3).coords(1,voxel), v.ROIs(3).coords(2,voxel)) = 1;
    end
    
    %get the voxels outside of brain and spread
    outsideMask = outsideMask - brainMask - spreadMask + 1;


%vector average and fit a cosine to the phase gradient
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
            y = ampMeanPhases(1:round(end/2));
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
            [maxr2, bestFrequencyIndex] = max(r2s);
            bestFrequency = testFrequencies(bestFrequencyIndex);
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

%do the overlay analyses - make histograms 
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
                
        
