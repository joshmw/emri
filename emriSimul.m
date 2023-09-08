function emriSimul(varargin)
% Purpose
%   Creates a ground-truth image time series, then creates new time series by sampling the original time series the way described in Toi et al. (line scans).
%   You can play around with adding and misaligning signals to the underlying movie to test how well the sampling approach works in theory.
% 
% Usage
%   emriSimul(varargin)
%
% Required
%    N/A
%
% Optional key/value pairs
%   'xdim': size of the image (x by x). This needs to be odd the way the simulation is set up to sample conjugate pairs.
%   'brainSize': size of the square brain (centered in image)
%   'signalType': Specifices type of systematic modulation. Right now, only supports harmonic (sine).
%   'addPinkNoise': 1 to add independent pink noise to every line acquisition.
%   'pinkNoiseStd': The standard deviation of the pink noise added.
%   'noiseLevel': add gaussian noise of that std to all inImage timepoints
%   'ntimePointsMs': length of one line acquisition
%   'sampleTimeMs': TR length. Takes a sample (1 line) every n miliseconds. i.e. nTimePointsMs = 500 and sampleTimeMS = 5 -> 100 samples of each line.
%   'hz': frequency of simulated signal in brain voxels (cycles/1000 ms)
%   'encodingDirection': the direction in wich you sample the individual lines. The specified direction is the FREQUENCY encoding direction, or the
%    fast acquisition direction.
%   'amplitude': amplitude of the cosine signal modulation in the brain voxels
%   'brainBaseContrast': base contrast of the brain (others pixels are 0).  
%   'phaseShiftBetweenLines': In the harmonic modulation case, this will randomize the phase of the modulation between each line acquisition.
%   'phaseShfitBetweenLinesMagnitude': Amount of randomization between line acquisitions. 1 = totally random, 0 = the same.
%   'pinkNoiseStd': 
%
% Output
%   N/A
%
% Description
%  
%  There are 4 main steps:
%   1. Create an "inImage", which is your ground-truth time series that you reconstruct from. This image will have a "brain", which signals go into. Other voxels are 0.
%   2. Mimic the Toi et al. procedure by going through the movie and, at each TR, fourier transform the image and grab single lines of the frequency
%   image at a time. The image is n cells (n=number of lines) of dimension xdim, ydim, ntimePointsMs. Each line is sampled sequentially from one cell.
%   This happens:
%       I. First, by using the first cell to sample frequency lines from. This is the same as assuming perfect stationarity and should give you the
%       original image (from cell 1). This creates "ksp".
%       II. Second, by using all the cells. You can randomize the phase of the signals in each cell to mimic non-stationarity. For example, some alpha
%       waves present during the entire acquisition will not be in perfect phase between each line collection. This creates "phaseUnlockedKsp".
%   3. Try and align the "phaseUnlockedKsp" from step 4. To do so, create movies from individual lines and try to align temporally in time by take the
%   absolute value of the fourier transform along the time dimension. There are other steps explained in the code. This creates "realignedKsp".
%   We don't really use this, but it's interesting to see. You can realign single component modulation to be cosine phase.
%   4. Plot results. As of now, the script will plot the image at a single timepoint using ksp (ground-truth), phaseUnlockedKsp, and realignedKsp. It will also plot
%   example voxel time series from the reconstructed movies. Finally, it will plot voxels outside of the "brain" along the phase-encoding direction.
%   
%   The flow of scripts goes: emriSimul -> inImageCreate -> inImageReconPhaseLocked -> inImageReconPhaseUnlocked -> realignUnlockedKsp.
%
%   Created by Josh Wilson some time in early 2023. Email joshmw@stanford.edu for questions.
%
%

%% Read input key/value pairs %%
p = inputParser;
%size of image and brain
p.addParameter('xdim',33,@isnumeric); % pixels x
p.addParameter('brainSize',8,@isnumeric);
%type of underlying noise signal. if harmonic, set the amplitude/frequency of modulation and line offset. if pink, set std.
p.addParameter('signalType','sine',@ischar)
p.addParameter('addPinkNoise',0,@islogical)
p.addParameter('amplitude',1,@isnumeric); %amplitude of hz modulation in sine case
p.addParameter('hz',3,@isnumeric); %signal frequency in "brain" voxels
p.addParameter('phaseShiftBetweenLines',1,@islogical); % will randomize the phase between lines instead of making continuous 
p.addParameter('phaseShiftBetweenLinesMagnitude',1,@isnumeric); %amount of shift between lines. 0 is none, 1 is fully random.
p.addParameter('pinkNoiseStd',0,@isnumeric); %amplitude of pink noise, if you are adding it
%set individual voxel offset if desired.
p.addParameter('numVoxelGroups',1,@isnumeric); %number of groups of voxels that oscillate together.
p.addParameter('offsetIndividualVoxels',1,@islogical); %1 to add a random offset to voxel phase in 'sine' condition
p.addParameter('offsetIndividualVoxelsMagnitude',1,@isnumeric); %how randomized individual voxel offset is. 0 is none (all same), 1 is fully random.
%add a visual response if desired. set the amount of shift, amplitude, and mean/std of peak.
p.addParameter('addSignal',0,@islogical); %add a gaussian signal
p.addParameter('signalShiftMagnitude',0,@isnumeric); %how much to shift the signal as % of the length of time series.
p.addParameter('signalAmp',.169,@isnumeric); % amplitude of the signal
p.addParameter('signalMean',100,@isnumeric);
p.addParameter('signalStd',15,@isnumeric);
%add gaussian noise to entire image if desired.
p.addParameter('noiseLevel',0,@isnumeric); %std of gaussian noise added to OG image
%other things you probably should not change - TR, encoding direction, length, luminance, etc.
p.addParameter('encodingDirection','x',@ischar); %FREQUENCY encoding direction (fast) - orthogonal to the phase-encoding direction
p.addParameter('ntimePointsMs',1000,@isnumeric); %time length of each line sample
p.addParameter('sampleTimeMs',5,@isnumeric); %TR length
p.addParameter('brainBaseContrast',100,@isnumeric); %Base brain value 
p.addParameter('buildConjugateLines',1,@isnumeric);  %Probably should not change. Half fourier approach.
p.addParameter('simNumber',1,@isnumeric);  %Probably should not change. Half fourier approach.
p.addParameter('graphStuff',1,@isnumeric);  %
p.addParameter('rng',100,@isnumeric); %rng seed - defaults to 100;

% make into parameters
p.parse;
p.addParameter('numKsamples',floor(p.Results.ntimePointsMs/p.Results.sampleTimeMs),@isnumeric)
p.parse(varargin{:}); params = p.Results;
params.numKsamples = floor(params.ntimePointsMs/params.sampleTimeMs);
params.ydim = params.xdim;
rng(params.rng) % set seed if you want to examine the effect of specific parameters

%% MAKE THE ORIGINAL IMAGE %%
inImage = inImageCreate(params);


%% GET THE DIFFERENT FREQUENCY IMAGES %%
%phase locked recon using only the first line sample of K space
ksp = inImageReconPhaseLocked(params, inImage);

% phase unlocked recon using all line samples
phaseUnlockedKsp = inImageReconPhaseUnlocked(params, inImage);

% try and align temporally...
realignedKsp = realignUnlockedKsp(params,phaseUnlockedKsp);

% get the true image sampled at the same frequency as kspace
trueImage = inImage{1}(:,:,1:params.sampleTimeMs:params.ntimePointsMs);

% reconstruct the images from k space
recoveredMovie = []; shiftRecoveredMovie = []; realignedRecoveredMovie = [];
for i = 1:params.numKsamples, recoveredMovie(:,:,i) = ifft2(ifftshift(ksp(:,:,i))); end
for i = 1:params.numKsamples, shiftRecoveredMovie(:,:,i) = ifft2(ifftshift(phaseUnlockedKsp(:,:,i))); end
for i = 1:params.numKsamples, realignedRecoveredMovie(:,:,i) = ifft2(ifftshift(realignedKsp(:,:,i))); end


%% PLOT THE RESULTS %%
if params.graphStuff
%pick an example voxel to look at, from the middle of the image
voxel = round(params.xdim/2);

%plot some example time series in the inImage
plotExampleVoxelTimeseries(params,inImage,voxel);

%pick a time point to compare static images
timePoint = round(rand*params.numKsamples);

figure; sgtitle('Example images and voxel time series')
subplot(3,3,1); imagesc(ifft2(ifftshift(ksp(:,:,timePoint)))); title('stationary sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar,
subplot(3,3,4); imagesc(abs(ksp(:,:,timePoint))); title('stationary sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)'),colorbar,

if ~params.buildConjugateLines
subplot(3,3,2); imagesc(abs(ifft2(ifftshift(phaseUnlockedKsp(:,:,timePoint))))); title('phase-shifted sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
else
subplot(3,3,2); imagesc(ifft2(ifftshift(phaseUnlockedKsp(:,:,timePoint)))); title('phase-shifted sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
end
subplot(3,3,5); imagesc(abs(phaseUnlockedKsp(:,:,timePoint))); title('phase-shifted sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)'),colorbar

subplot(3,3,3); imagesc(ifft2(ifftshift(realignedKsp(:,:,timePoint)))); title('phase-realigned sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
subplot(3,3,6); imagesc(abs(realignedKsp(:,:,timePoint))); title('phase-realigned sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)'),colorbar

subplot(3,3,7), title('Example brain voxel time series'), xlabel('Time'), ylabel('Magnitude'), hold on,
plot(squeeze(recoveredMovie(voxel,voxel,:)))
plot(squeeze(shiftRecoveredMovie(voxel,voxel,:)))
plot(squeeze(realignedRecoveredMovie(voxel,voxel,:)))
legend('Resampled line (ground truth)','Different line samples','Phase-realigned, different line samples')

subplot(3,3,8), hold on, title('Non-brain voxel timeseries (shifted lines)'), xlabel('Time'), ylabel('Magnitude'),
plotNonBrainVoxels(params,shiftRecoveredMovie);

subplot(3,3,9), hold on, title('Non-brain voxel timeseries (realigned from shifted lines)'), xlabel('Time'), ylabel('Magnitude'),
plotNonBrainVoxels(params,realignedRecoveredMovie);

end
keyboard
%%%%%% END OF SIMULATION %%%%%%%%

%save things if you want
%ts = reshape(shiftRecoveredMovie,[params.xdim params.ydim 1 params.numKsamples]); 
%[d h] = mlrImageLoad('testName.img');   
%str = sprintf('TR%2.0f',params.sampleTimeMs);
%str = replace(str, '.', '_');
%mlrImageSave(str,ts,h);
%close all

%saving on luxardo
ts = reshape(shiftRecoveredMovie,[params.xdim params.ydim 1 params.numKsamples]);
j = sprintf('numGroups1dim%4.4g',params.numVoxelGroups)
j = replace(j, '.', '_');
j = replace(j, ' ', '')
save(j,'ts','params');
close all

ts = reshape(realignedRecoveredMovie,[params.xdim params.ydim 1 params.numKsamples]);
save(strcat('realigned',j),'ts','params');
%for TR = [2 4 5 8 10 20 25 40 50 100], emriSimul('sampleTimeMs',TR),end





%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%% plotNonBrainVoxels %%
function plotNonBrainVoxels(params,image)

brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);

switch params.encodingDirection
    %x phase encoding direction
    case 'x'
        for x = 1:params.xdim;
            for y = 1:params.ydim;
                if (~ismember(x,brain) & ismember(y,brain));
                    plot(1:params.numKsamples,reshape(image(x,y,:),1,params.numKsamples))
                end
            end
        end
    
    %y phase encoding direction
    case 'y'
        for x = 1:params.xdim;
            for y = 1:params.ydim;
                if (ismember(x,brain) & ~ismember(y,brain));
                    plot(1:params.numKsamples,reshape(image(x,y,:),1,params.numKsamples))
                end
            end
        end
end
%ylim([-params.amplitude params.amplitude])


%% plotBrainVoxels %%
function plotBrainVoxels(params,image)

brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2);

    %x phase encoding direction
    for x = 1:params.xdim;
        for y = 1:params.ydim;
            if (ismember(x,brain) & ismember(y,brain));
                    plot(1:params.numKsamples,reshape(image(x,y,:),1,params.numKsamples))
            end
        end
    end




%% plotExampleVoxelTimeseries %%
function plotExampleVoxelTimeseries(params,inImage,voxel)
%for 2 example voxels, concatenate the time series from the first 3 cells
wholeSeries = []; wholeSeries2 = [];
for i  = 1:3, wholeSeries = [wholeSeries squeeze(inImage{i}(voxel,voxel,:))']; end     
for i  = 1:3, wholeSeries2 = [wholeSeries2 squeeze(inImage{i}(voxel+2,voxel+2,:))']; end

%plot the time series
figure, hold on,
plot(wholeSeries),plot(wholeSeries2);

%plot the breaks between line acquisitions
voxelMin = min(min(wholeSeries),min(wholeSeries2));
voxelMax = max(max(wholeSeries),max(wholeSeries2));
plot([params.ntimePointsMs params.ntimePointsMs], [voxelMin voxelMax],'lineWidth',2,'color','k')
plot([params.ntimePointsMs*2 params.ntimePointsMs*2], [voxelMin voxelMax],'lineWidth',2,'color','k')
xlabel('Time (ms)'),ylabel('Voxel magnitude'),legend('Example voxel 1','Example voxel 2','',''),title('Example voxel time series in the inImage')

%label the figure
xticks([0 params.ntimePointsMs/2 params.ntimePointsMs params.ntimePointsMs*1.5 params.ntimePointsMs*2 params.ntimePointsMs*2.5 params.ntimePointsMs*3]) 
xticklabels({0, 'Line 1 acquistition time series', params.ntimePointsMs, 'Line 2 acquisition time series', params.ntimePointsMs*2, 'Line 3 acquisition time series', params.ntimePointsMs*3}) 












function extraStuff

% saving a simulation
ts = reshape(shiftRecoveredMovie,[params.xdim params.ydim 1 params.numKsamples]); 
[d h] = mlrImageLoad('testName.img');   
mlrImageSave('33p3hzHomogenousLineOffsetYEncode',ts,h);


%running mrClassify with mrsession open
mrClassify(params)



%old
figure, hold on
    for col = brain
        plot(1:params.ydim,cos(v.analyses{1}.overlays(frequencyToView*3).data{scan}(:,col)))
        plot(1:params.ydim,sin(v.analyses{1}.overlays(frequencyToView*3).data{scan}(:,col)),'--')

    end

figure, hold on
    for row = 1:params.xdim
        plot(1:length(brain),cos(v.analyses{1}.overlays(frequencyToView*3).data{scan}(row,brain)))
        plot(1:length(brain),sin(v.analyses{1}.overlays(frequencyToView*3).data{scan}(row,brain)),'--')
    end
