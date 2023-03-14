%% emriSimul.m %%
%
% Simluation for a spectral emri analysis. Synthesize some "brain" voxels
% with specific modulation frequencies, then image those with the emri
% sequence (take the same line of k space many times, then take the next
% line many times, etc). 
%
% The original image is generated as an array of matrices, where each
% element of the array has the XYT time series for each line scan.
% Voxel{n+1}(x,y,t=1) will pick up where Voxel{n}(x,y,t=end of line scan)
% leaves off, phase wise, so concat( voxel{n}(x,y,:), voxel{n+1}(x,y,:) )
% should be continuous. If the signal frequency is uneven with the length
% of each line sample, you get signal dephasing between lines. 
% 
% The "True image" is just the first {n} of the original image. The reconstructed
% True image is reconstructed from {n},{n},...{n} rather than
% {n},{n+1},...{n end} - this is the "shifted recon" image, akin to what we actually measure. 
% 
% I am working on a phase-aligned reconstruction, but it is not ready yet.
% I will make comments on code where to ignore that.
%
% PARAMETERS:
%   'xdim': number of voxels in the whole image
%   'brainSize': size of the brain (centered in image)
%   'noiseLevel': add gaussian noise of that std to all OG image timepoints
%   'ntimePointsMs': length of one line acquisition.
%   'hz': frequency of simulated signal in brain voxels (cycles/1000 ms)
%   'encodingDirection': the direction in wich you sample the individual lines
%   'amplitude': amplitude of the cosine signal modulation in the brain voxels
%   'sampleTimeMs': TR length. Takes a sample (1 line) every n miliseconds. i.e. nTimePointsMs = 500 and sampleTimeMS = 5 -> 100 samples
%   'brainBaseContrast': base contrast of the brain (others are 0).  


function emriSimul(varargin)

%% set parameters %%
xdim = 32; %pixels x
brainSize = 16;
noiseLevel = 0; %std of gaussian noise added to OG image
ntimePointsMs = 500 %length of simulation
encodingDirection = 'x'; %phase encoding direction
hz = 3.99; %signal frequency in "brain" voxels
amplitude = 1;
sampleTimeMs = 5; %TR length
brainBaseContrast = 0;
numKsamples = ntimePointsMs/sampleTimeMs; %number of samples will be total length/TR length

getArgs(varargin);
ydim = xdim; %pixels y


%% make an image with modulatory signals %%
originalImage = makeOriginalImage(xdim,ydim,noiseLevel,ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz)


%% get the different reconstructions
%phase locked recon using only the first line sample of K space
ksp = doPhaseLockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);

% phase unlocked recon using ann line samples
phaseUnlockedKsp = doPhaseUnlockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);

% try and phase correct the lines in the unlocked condition - THIS DOES NOT WORK YET!!
realignedKsp = doCorrectedPhaseRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);


%% plot things %%
% get the true image sampled at the same frequency as kspace
trueImage = originalImage{1}(:,:,1:sampleTimeMs:ntimePointsMs);

%pick a time point to compare
timePoint = round(rand*ntimePointsMs/sampleTimeMs);

figure;
subplot(2,4,1); imshow(trueImage(:,:,timePoint)); title('ground truth image'); xlabel('x (image)'); ylabel('y (image)')
subplot(2,4,5); imshow(fft2(trueImage(:,:,timePoint))); title('ground truth kspace'); xlabel('x (frequency)'); ylabel('y (frequency)')

subplot(2,4,2); imshow(ifft2(ksp(:,:,timePoint))); title('stationary sampled image'); xlabel('x (image)'); ylabel('y (image)')
subplot(2,4,6); imshow(ksp(:,:,timePoint)); title('stationary sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)')

subplot(2,4,3); imshow(ifft2(phaseUnlockedKsp(:,:,timePoint))); title('phase-shifted sampled image'); xlabel('x (image)'); ylabel('y (image)')
subplot(2,4,7); imshow(phaseUnlockedKsp(:,:,timePoint)); title('phase-shifted sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)')

%subplot(2,4,4); imshow(ifft2(realignedKsp(:,:,timePoint))); title('phase-realigned sampled image - DOESNT WORK RIGHT NOW'); xlabel('x (image)'); ylabel('y (image)')
%subplot(2,4,8); imshow(realignedKsp(:,:,timePoint)); title('phase-realigned sampled kspace - DOESNT WORK RIGHT NOW'); xlabel('x (frequency)'); ylabel('y (frequency)')

% reconstruct the images from k space
recoveredImage = [];
for i = 1:numKsamples,
    recoveredImage(:,:,i) = ifft2(ksp(:,:,i));
end

shiftRecoveredImage = [];
for i = 1:numKsamples,
    shiftRecoveredImage(:,:,i) = ifft2(phaseUnlockedKsp(:,:,i));
end

realignedRecoveredImage = [];
for i = 1:numKsamples,
    realignedRecoveredImage(:,:,i) = ifft2(realignedKsp(:,:,i));
end

% plot example voxels: an original image one, a phase-alligned recon one,
% a phase-shifted recon one, and a "false positive" from the phase shifted recon
h = trueImage(xdim/2,ydim/2,:);
j = recoveredImage(xdim/2,ydim/2,:);
k = shiftRecoveredImage(xdim/2,ydim/2,:);
l = shiftRecoveredImage(1,ydim/2,:);;

figure,title('Example voxels')
subplot(2,2,1);plot(1:sampleTimeMs:ntimePointsMs,h(:)); title('Ground truth'); ylim([-1 1]); xlabel('time (ms)'); ylabel('signal (a.u.)');
subplot(2,2,2);plot(1:sampleTimeMs:ntimePointsMs,j(:)); title('Phase-locked measurement recon'); ylim([-1 1]); xlabel('time (ms)'); ylabel('signal (a.u.)');
subplot(2,2,3);plot(1:sampleTimeMs:ntimePointsMs,k(:)); title('Phase-scrambled measurement recon'); ylim([-1 1]); xlabel('time (ms)'); ylabel('signal (a.u.)');
subplot(2,2,4);plot(1:sampleTimeMs:ntimePointsMs,l(:)); title('Phase-scrambled recon (outside brain)'); ylim([-1 1]); xlabel('time (ms)'); ylabel('signal (a.u.)');

sgtitle('Example voxels')

%%%% END OF SIMULATION %%%%














%%%%%%%%%%%%%%%%%%%%%%%
%% makeOriginalImage %%
%%%%%%%%%%%%%%%%%%%%%%%
function originalImage = makeOriginalImage(xdim,ydim,noiseLevel,ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz)

%make the image with a "brain" which will have signals in it
for lineSample = 1:xdim
    originalImage{lineSample} = normrnd(0,noiseLevel,xdim,ydim,ntimePointsMs);
end

%put a signal of determined frequency but random phase in the "brain" voxels
brain = (round(xdim/2)-brainSize/2):(round(xdim/2)+brainSize/2); %brain is square here, and xdim=ydim

for voxelx = 1:xdim;
    for voxely = 1:ydim;
        if ismember(voxelx,brain) & ismember(voxely,brain)
            voxelOffset = round(rand*ntimePointsMs); %makes the voxels different phase from eachother
            %voxelOffset = 0;
            parfor lineSample = 1:xdim
                originalImage{lineSample}(voxelx,voxely,:) = brainBaseContrast + amplitude*cos( ((1+(lineSample-1)*ntimePointsMs):(lineSample*ntimePointsMs))*2*pi/1000*hz+voxelOffset);
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%
%% doNormalRecon %%
%%%%%%%%%%%%%%%%%%%

function ksp = doPhaseLockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection)

%initialize kspace for speed
ksp = zeros(xdim,ydim,numKsamples);

% sample each line (row or column) at each sample timepoint (TR length)
for row = 1:xdim
    %for each TR, sample the line in kspace and add to the kspace image
    for sample = 1:numKsamples  
        imageKsp = fft2(originalImage{1}(:,:,sampleTimeMs*(1+(sample-1))));      
        switch encodingDirection
            case 'x'
                kspLine = imageKsp(row,:);
                ksp(row,:,sample) = kspLine;
            case 'y'
                kspLine = imageKsp(:,row);
                ksp(:,row,sample) = kspLine;
        end
    end
end



%%%%%%%%%%%%%%%%%%%%
%% doShiftedRecon %%
%%%%%%%%%%%%%%%%%%%%

% do the same as above, but shifting the phase of the voxels
function phaseUnlockedKsp = doPhaseUnlockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);

for row = 1:xdim
    for sample = 1:numKsamples
        %takes the temporal order of the images
        imageKsp = fft2(originalImage{row}(:,:,sampleTimeMs*(1+(sample-1))));      
        switch encodingDirection
            case 'x'
                kspLine = imageKsp(row,:);
                phaseUnlockedKsp(row,:,sample) = kspLine;
            case 'y'
                kspLine = imageKsp(:,row);
                phaseUnlockedKsp(:,row,sample) = kspLine;   
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% doCorrectedPhaseRecon %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function realignedKsp = doCorrectedPhaseRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);

unalignedKsp = zeros(xdim,xdim,numKsamples);
realignedKsp = zeros(xdim,xdim,numKsamples);

% get the kspace image the way they did originally
for row = 1:xdim
    %shift the signals a random amount
    for sample = 1:numKsamples
        imageKsp = fft2(originalImage{row}(:,:,sampleTimeMs*(1+(sample-1))));      
        switch encodingDirection
            case 'x'
                kspLine = imageKsp(row,:);
                unalignedKsp(row,:,sample) = kspLine;
                
            case 'y'
                kspLine = imageKsp(:,row);
                unalignedKsp(:,row,sample) = kspLine;
        end
    end
    realignedLine = phaseAlignKsp(unalignedKsp,row,encodingDirection,numKsamples);
    realignedKsp(row,:,:) = realignedLine;
end


%%%%%%%%%%%%%%%%%%%
%% phaseAlignKsp %%
%%%%%%%%%%%%%%%%%%%
function realignedLine = phaseAlignKsp(unalignedKsp,row,encodingDirection,numKsamples)

%get just the single line
singleLineFxFyT = zeros(size(unalignedKsp));
singleLineFxFyT(row,:,:) = unalignedKsp(row,:,:);

%recon an XYT image from the FxFyT matrix
for sample = 1:numKsamples,
    singleLineXYT(:,:,sample) = ifft2(singleLineFxFyT(:,:,sample));
end

%take the fourier transform in the time dimension
singleLineXYFt = fft(singleLineXYT,[],3);

%take the absolute value to throw out temporal phase information
%absSingleLineXYFt = abs(singleLineXYFt);
absSingleLineXYFt = singleLineXYFt;

%go back to XYT with temporal phase thrown out
nophaseSingleLineXYT = ifft(absSingleLineXYFt,[],3);

%go back to FxFyT now that you've removed temporal phase
for sample = 1:numKsamples,
    nophaseSingleLineFxFyT(:,:,sample) = fft2(nophaseSingleLineXYT(:,:,sample));
end
realignedLine = nophaseSingleLineFxFyT(row,:,:);


