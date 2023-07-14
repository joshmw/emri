function emriGradientAnalPlot(SIDs)
%calling function for emriGradientAnal. Will cycle through all subjects and plot the stuff you want.
% emriGradientAnalPlot(["s0615" "s0616" "s0617" "s0618" "s0619"])

%go to data directory
cd ~/data/EMRI/paperData/

%% get the data you need

% loop over scans. 1 = gated, 2 = ungated,, 3 = flipped.
for scan = 1
    %loop over subjects
    for subject = 1:length(SIDs)
        % we dont have flipped data for subject s0615.
        if ~((SIDs(subject) == 's0615') & scan==3)
            cd(SIDs(subject));
            %load the view, ROIs, Analyses, and go to ungated scan
            v = newView;
            v = viewSet(v,'curGroup','Averages')
            v = viewSet(v,'curScan',scan)
            v = loadAnalysis(v,'emriAnal');
            v = loadROI(v,'spread') ;
            v = loadROI(v,'insideBrain') ;
            v = loadROI(v,'outsideBrain') ;
            if ~(SIDs(subject) == 's0615')
                v = loadROI(v,'spreadFlipped') ;
                v = loadROI(v,'insideBrainFlipped') ;
                v = loadROI(v,'outsideBrainFlipped') ;
            end
            %get the frequencies and amplitudes from emriGradientAnal
            [brainMean, spreadMean, outsideMean, brainStd, spreadStd, outsideStd, brainNum, spreadNum, outsideNum] = emriGradientAnal(v);  
            %save the data into a larger data table
             data(subject,:) = [brainMean, spreadMean, outsideMean, brainStd, spreadStd, outsideStd, brainNum, spreadNum, outsideNum];
            %clear the view before going back and doing for each subject
            pause(1)
            deleteView(v)
            mrQuit
            cd ..
        else
            %if we are on subject 1, scan 3
            data(subject,:) = NaN(1,9);
        end

    end
   
    %% once we get all the data for the scan type, graph.

    figure(150), subplot(1,3,scan), hold on
    %loop over subjects
    for sub = 1:height(data)
        scatter(sub,data(sub,1),100,'MarkerFaceColor','g','MarkerEdgeColor','w')
        scatter(sub+20,data(sub,2),100,'MarkerFaceColor','c','MarkerEdgeColor','w')
        scatter(sub+40,data(sub,3),100,'MarkerFaceColor','m','MarkerEdgeColor','w')
        %errorbars, if you want
        %errorbar(sub,data(sub,1),data(sub,4)/data(sub,7),'Color','g')
        %errorbar(sub+20,data(sub,2),data(sub,5)/data(sub,8),'Color','c')
        %errorbar(sub+40,data(sub,3),data(sub,6)/data(sub,9),'Color','m')
        ylim([0 1]),xlim([0 50])
    end

    %plot the group average autocor inside the brain
    scatter(median(1:height(data)), mean(data(:,1), 'omitnan'), 150, 'filled', 'k')
    errorbar(median(1:height(data)), mean(data(:,1), 'omitnan'), std(data(:,1), 'omitnan')/height(data),'k','MarkerEdgeColor','w'), 
    %plot the group average autocor along the spread dimension
    scatter(20+median(1:height(data)), mean(data(:,2), 'omitnan'), 150, 'filled', 'k','MarkerEdgeColor','w')
    errorbar(20+median(1:height(data)), mean(data(:,2), 'omitnan'), std(data(:,2), 'omitnan')/height(data),'k'), 
    %plot the group average autocor outside the spread
    scatter(40+median(1:height(data)), mean(data(:,3), 'omitnan'), 150, 'filled', 'k')
    errorbar(40+median(1:height(data)), mean(data(:,3), 'omitnan'), std(data(:,3), 'omitnan')/height(data),'k','MarkerEdgeColor','w'), 
    
end

%label each subplot with what kind of data it has
subplot(1,3,1), title('Cardiac Gated')
subplot(1,3,2), title('Ungated')
subplot(1,3,3), title('flipped')
keyboard

