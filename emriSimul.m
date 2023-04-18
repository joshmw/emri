function emriSimul(varargin)
% One line description of the function here
% 
% Synopsis
%   emriSimuli(varargin)
%
% Required
%    N/A
%
% Optional key/value pairs
%   'xdim': number of voxels in the whole image
%   'brainSize': size of the brain (centered in image)
%   'noiseLevel': add gaussian noise of that std to all OG image timepoints
%   'ntimePointsMs': length of one line acquisition.
%   'hz': frequency of simulated signal in brain voxels (cycles/1000 ms)
%   'encodingDirection': the direction in wich you sample the individual lines
%   'amplitude': amplitude of the cosine signal modulation in the brain voxels
%   'sampleTimeMs': TR length. Takes a sample (1 line) every n miliseconds. i.e. nTimePointsMs = 500 and sampleTimeMS = 5 -> 100 samples
%   'brainBaseContrast': base contrast of the brain (others are 0).  
%
% Output
%  N/A
%
% Description
%
%  Simluation for a spectral emri analysis. Synthesize some "brain"
%  voxels with specific modulation frequencies, then image those with
%  the emri sequence (take the same line of k space many times, then
%  take the next line many times, etc).
%
%  The original image is generated as an array of matrices, where each
%  element of the array has the XYT time series for each line scan.
%  Voxel{n+1}(x,y,t=1) will pick up where Voxel{n}(x,y,t=end of line
%  scan) leaves off, phase wise, so concat( voxel{n}(x,y,:),
%  voxel{n+1}(x,y,:) ) should be continuous. If the signal frequency
%  is uneven with the length of each line sample, you get signal
%  dephasing between lines.
% 
% The "True image" is just the first {n} of the original image. The
% reconstructed True image is reconstructed from {n},{n},...{n} rather
% than {n},{n+1},...{n end} - this is the "shifted recon" image, akin
% to what we actually measure.
% 
% I am working on a phase-aligned reconstruction, but it is not ready yet.
% I will make comments on code where to ignore that.
%
% See also
%


%% Read input key/value pairs

p = inputParser;
p.addParameter('xdim',17,@isnumeric);
p.addParameter('ydim',17,@isnumeric);
p.addParameter('brainSize',4,@isnumeric);
p.addParameter('encodingDirection','x',@ischar);
p.addParameter('offsetIndividualVoxels',0,@islogical);
p.addParameter('phaseShiftBetweenLines',1,@islogical);
p.addParameter('noiseLevel',0,@isnumeric);
p.addParameter('ntimePointsMs',1000,@isinteger);
p.addParameter('hz',4,@isnumeric);
p.addParameter('amplitude',1,@isnumeric);
p.addParameter('sampleTimeMs',5,@isnumeric);
p.addParameter('brainBaseContrast',1000,@isnumeric);
p.addParameter('buildConjugateLines',1,@isnumeric);  % Probably should not change
p.parse;

% Build the parameters

xdim = p.Results.xdim; % pixels x
ydim = p.Results.ydim; % pixels y
brainSize         = p.Results.brainSize;
encodingDirection = p.Results.encodingDirection; %phase encoding direction
offsetIndividualVoxels = p.Results.offsetIndividualVoxels; %1 to add a random offset to voxel phase in 'sine' condition
phaseShiftBetweenLines = p.Results.phaseShiftBetweenLines; % will randomize the phase between lines instead of making continuous 
noiseLevel = p.Results.noiseLevel;       %std of gaussian noise added to OG image
ntimePointsMs = p.Results.ntimePointsMs; %length of simulation
hz = p.Results.hz;                       %signal frequency in "brain" voxels
amplitude = p.Results.amplitude;
sampleTimeMs = p.Results.sampleTimeMs; %TR length
brainBaseContrast = p.Results.brainBaseContrast;
buildConjugateLines = p.Results.buildConjugateLines;

numKsamples = ntimePointsMs/sampleTimeMs; %number of samples will be total length/TR length

% JLG Group method
% getArgs(varargin);

%% MAKE THE ORIGINAL IMAGE %%
%
%  nounVerb 
%   inImage.create
%   inImage.set
%   inImage.get
%   inImage.resample
%
originalImage = makeOriginalImageSine(xdim,ydim,noiseLevel,ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz,offsetIndividualVoxels,phaseShiftBetweenLines);


%% GET THE DIFFERENT RECONSTRUCTINONS %%
%phase locked recon using only the first line sample of K space
ksp = doPhaseLockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);

% phase unlocked recon using ann line samples
phaseUnlockedKsp = doPhaseUnlockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection,buildConjugateLines);

% try and align temporally...
realignedKsp = correctPhaseUnlockedKsp(phaseUnlockedKsp,xdim,numKsamples,encodingDirection);


%% PLOT THE RESULTS %%
% get the true image sampled at the same frequency as kspace
trueImage = originalImage{1}(:,:,1:sampleTimeMs:ntimePointsMs);

% reconstruct the images from k space
recoveredImage = []; shiftRecoveredImage = []; realignedRecoveredImage = [];
for i = 1:numKsamples, recoveredImage(:,:,i) = ifft2(ifftshift(ksp(:,:,i))); end
for i = 1:numKsamples, shiftRecoveredImage(:,:,i) = ifft2(ifftshift(phaseUnlockedKsp(:,:,i))); end
for i = 1:numKsamples, realignedRecoveredImage(:,:,i) = ifft2(ifftshift(realignedKsp(:,:,i))); end

%pick a time point to compare
timePoint = round(rand*ntimePointsMs/sampleTimeMs)

figure; sgtitle('Example images and voxel time series')
%subplot(3,4,1); imagesc(trueImage(:,:,timePoint)); title('ground truth image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
%subplot(3,4,5); imagesc(abs(fftshift(fft2(trueImage(:,:,timePoint))))); title('ground truth kspace'); xlabel('x (frequency)'); ylabel('y (frequency)') 

subplot(3,3,1); imagesc(ifft2(ifftshift(ksp(:,:,timePoint)))); title('stationary sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
subplot(3,3,4); imagesc(abs(ksp(:,:,timePoint))); title('stationary sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)')

if ~buildConjugateLines
subplot(3,3,2); imagesc(abs(ifft2(ifftshift(phaseUnlockedKsp(:,:,timePoint))))); title('phase-shifted sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
else
subplot(3,3,2); imagesc(ifft2(ifftshift(phaseUnlockedKsp(:,:,timePoint)))); title('phase-shifted sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
end
subplot(3,3,5); imagesc(abs(phaseUnlockedKsp(:,:,timePoint))); title('phase-shifted sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)')

subplot(3,3,3); imagesc(ifft2(ifftshift(realignedKsp(:,:,timePoint)))); title('phase-realigned sampled image'); xlabel('x (image)'); ylabel('y (image)'), colormap(gray), colorbar
subplot(3,3,6); imagesc(abs(realignedKsp(:,:,timePoint))); title('phase-realigned sampled kspace'); xlabel('x (frequency)'); ylabel('y (frequency)')

subplot(3,3,7), title('Example brain voxel time series'), xlabel('Time'), ylabel('Magnitude'), hold on, voxel = round(xdim/2);
plot(squeeze(recoveredImage(voxel,voxel,:)))
plot(squeeze(shiftRecoveredImage(voxel,voxel,:)))
plot(squeeze(realignedRecoveredImage(voxel,voxel,:)))
legend('Resampled line (ground truth)','Different line samples','Phase-realigned, different line samples')

subplot(3,3,8), hold on, title('Non-brain voxel timeseries (shifted lines)'), xlabel('Time'), ylabel('Magnitude'),
plotNonBrainVoxels(xdim,ydim,numKsamples,shiftRecoveredImage,brainSize,encodingDirection);

subplot(3,3,9), hold on, title('Non-brain voxel timeseries (realigned from shifted lines)'), xlabel('Time'), ylabel('Magnitude'),
plotNonBrainVoxels(xdim,ydim,numKsamples,realignedRecoveredImage,brainSize,encodingDirection);


%%%%%% END OF SIMULATION %%%%%%%%
keyboard

% sz = size(shiftRecoveredImage);
% niftiImage = reshape(abs(shiftRecoveredImage),[sz(1) sz(2) 1 sz(3)]);











%% HELPER FUNCTIONS %%

%%%%%%%%%%%%%%%%%%%%%%%
%% makeOriginalImageSine %%
%%%%%%%%%%%%%%%%%%%%%%%
function originalImage = makeOriginalImageSine(xdim,ydim,noiseLevel,ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz,offsetVoxels,phaseShiftBetweenLines);
%
% 
% Make the image with a "brain" which will have signals in it
originalImage = cell(xdim,1);
for lineSample = 1:xdim
    originalImage{lineSample} = normrnd(0,noiseLevel,xdim,ydim,ntimePointsMs);
end

%calculate line offsets if randomizing
if phaseShiftBetweenLines 
    lineOffsetMs = round(rand(1,xdim)*ntimePointsMs);
end

%put a signal of determined frequency but random phase in the "brain" voxels
brain = (round(xdim/2)-brainSize/2):(round(xdim/2)+brainSize/2); %brain is square here, and xdim=ydim

for voxelx = 1:xdim
    for voxely = 1:ydim
        if ismember(voxelx,brain) && ismember(voxely,brain)

            %offset voxel phase if desired
            if offsetVoxels
                voxelOffset = round(rand*ntimePointsMs); %makes the voxels different phase from eachother
            else
                voxelOffset = 0;
            end

            %make the lines
            if ~phaseShiftBetweenLines
                parfor lineSample = 1:xdim
                    originalImage{lineSample}(voxelx,voxely,:) = brainBaseContrast + amplitude*cos( ((1+(lineSample-1)*ntimePointsMs):(lineSample*ntimePointsMs))*2*pi/1000*hz+voxelOffset);
                end    
            else
                parfor lineSample = 1:xdim
                    originalImage{lineSample}(voxelx,voxely,:) = brainBaseContrast + amplitude*cos( ((1+lineOffsetMs(lineSample)):(lineOffsetMs(lineSample))+ntimePointsMs)*2*pi/1000*hz+voxelOffset);
                end
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
        imageKsp = fftshift(fft2(originalImage{1}(:,:,1+sampleTimeMs*(sample-1))));      
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
function phaseUnlockedKsp = doPhaseUnlockedRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection,buildConjugateLines);

%can sample 2 lines at once if you'd like to make kspace conjugate symmetric.
%this ensures a REAL IMAGE. effectively the same as half-fourier sampling.
if buildConjugateLines

    for row = 1:(floor(xdim/2)+1)
        for sample = 1:numKsamples
            %takes the temporal order of the images
            imageKsp = fftshift(fft2(originalImage{row}(:,:,1+sampleTimeMs*(sample-1))));      
            switch encodingDirection
                case 'x'
                    kspLine = imageKsp(row,:);
                    phaseUnlockedKsp(row,:,sample) = kspLine;
                    kspLine2 = imageKsp(xdim-row+1,:);
                    phaseUnlockedKsp(xdim-row+1,:,sample) = kspLine2;
                case 'y'
                    kspLine = imageKsp(:,row);
                    phaseUnlockedKsp(:,row,sample) = kspLine;
                    kspLine2 = imageKsp(:,xdim-row+1);
                    phaseUnlockedKsp(:,xdim-row+1,sample) = kspLine2;
            end
        end
    end
    
else
    %sampling one line at a time, no conjugate pairing.
    for row = 1:(xdim)
        for sample = 1:numKsamples
            %takes the temporal order of the images
            imageKsp = fftshift(fft2(originalImage{row}(:,:,1+sampleTimeMs*(sample-1))));      
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

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correctPhaseUnlockedKsp %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function realignedKsp = correctPhaseUnlockedKsp(phaseUnlockedKsp,xdim,numKsamples,encodingDirection)

realignedKsp = zeros(size(phaseUnlockedKsp));

for row = 1:(floor(xdim/2)+1)

    %take a single pair of FxFyT lines
    singleLineFxFyT = zeros(size(phaseUnlockedKsp));
    switch encodingDirection
        case 'x'
            singleLineFxFyT(row,:,:) = phaseUnlockedKsp(row,:,:);
            singleLineFxFyT(xdim-row+1,:,:) = phaseUnlockedKsp(xdim-row+1,:,:);
        case 'y'
            singleLineFxFyT(:,row,:) = phaseUnlockedKsp(:,row,:);
            singleLineFxFyT(:,xdim-row+1,:) = phaseUnlockedKsp(:,xdim-row+1,:);
    end

    %recon an XYT image from the single line FxFyT matrix
    for sample = 1:numKsamples,
    singleLineXYT(:,:,sample) = ifft2(ifftshift(singleLineFxFyT(:,:,sample)));
    end

    %take the fourier transform in the time dimension
    singleLineXYFt = fft(singleLineXYT,[],3);

    %take the absolute value to throw out temporal phase information, then flip back to original sign w/ DC matrix
    DCmatrix = (singleLineXYFt(:,:,1)>0)*2-1;
    absSingleLineXYFt = abs(singleLineXYFt);
    absSingleLineXYFtDC = absSingleLineXYFt .* DCmatrix;

    %go back to XYT with temporal phase thrown out
    nophaseSingleLineXYT = ifft(absSingleLineXYFtDC,[],3);

    %go back to FxFyT now that you've removed temporal phase
    for sample = 1:numKsamples,
        nophaseSingleLineFxFyT(:,:,sample) = fftshift(fft2(nophaseSingleLineXYT(:,:,sample)));
    end
    
    for sample = 1:numKsamples
        realignedKsp(:,:,sample) = realignedKsp(:,:,sample) + nophaseSingleLineFxFyT(:,:,sample);
        %realignedKsp(row,:,sample) = nophaseSingleLineFxFyT(row,:,sample);
        %realignedKsp(xdim-row+1,:,sample) = nophaseSingleLineFxFyT(xdim-row+1,:,sample);
    end
  
    if row == floor(xdim/2);
        plotPhaseAlignment(row,xdim,singleLineFxFyT,singleLineXYT,nophaseSingleLineFxFyT,nophaseSingleLineXYT);
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%
%% plotPhaseAlignment %%
%%%%%%%%%%%%%%%%%%%%%%%%
%this only works for 'x' phase encoding direction - could make it work
%otherwise, but I am lazy. will add later.
function plotPhaseAlignment(row,xdim,singleLineFxFyT,singleLineXYT,nophaseSingleLineFxFyT,nophaseSingleLineXYT)

fig = mlrSmartfig('line image')

%plot a snapshot of the frequency image
subplot(2,4,1)
imagesc(abs(singleLineFxFyT(:,:,1))),colorbar,colormap(gray)
xlabel('Fx'),ylabel('Fy'),title('OG Frequency image (1 timepoint)')

%plot the timecourse of the individual frequencies
subplot(4,4,2), hold on,
for x = 1:xdim
    plot(1:200,reshape(real(singleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Real Value Frequency component magnitudes over time (1 row)')

subplot(4,4,6), hold on,
for x = 1:xdim
    plot(1:200,reshape(imag(singleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Imaginary Value Frequency component magnitudes over time (1 row)')

%plot a snapshot of the actual image
subplot(2,4,3);
imagesc(singleLineXYT(:,:,1));colorbar,colormap(gray)
xlabel('x'),ylabel('y'),title('image created from single lines')

%plot the timecourse of the image pixels
subplot(2,4,4), hold on
for x = 1:xdim
plot(1:200,reshape(singleLineXYT(x,8,:),1,200))
end
xlabel('time'),ylabel('magnitude'),title('pixel values over time (1 column)')

%%then do all of that for the phase-corrected image...
subplot(2,4,5)
imagesc(abs(nophaseSingleLineFxFyT(:,:,1))),colorbar,colormap(gray)
xlabel('Fx'),ylabel('Fy'),title('OG Frequency image (1 timepoint)')

subplot(4,4,10), hold on,
for x = 1:xdim
plot(1:200,reshape(real(nophaseSingleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Real Value Frequency component magnitudes over time (1 row)')

subplot(4,4,14), hold on,
for x = 1:xdim
plot(1:200,reshape(imag(nophaseSingleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Imaginary Value Frequency component magnitudes over time (1 row)')

subplot(2,4,7);
imagesc(nophaseSingleLineXYT(:,:,1));colorbar,colormap(gray)
xlabel('x'),ylabel('y'),title('Realigned image')

subplot(2,4,8), hold on
for x = 1:xdim
plot(1:200,reshape(nophaseSingleLineXYT(x,8,:),1,200))
end
xlabel('time'),ylabel('magnitude')



%%%%%%%%%%%%%%%%%%%%%%%%
%% plotNonBrainVoxels %%
%%%%%%%%%%%%%%%%%%%%%%%%
function plotNonBrainVoxels(xdim,ydim,numKsamples,data,brainSize,encodingDirection)

brain = (round(xdim/2)-brainSize/2):(round(xdim/2)+brainSize/2);

switch encodingDirection
    %x phase encoding direction
    case 'x'
        for x = 1:xdim;
            for y = 1:ydim;
                if (~ismember(x,brain) & ismember(y,brain));
                    plot(1:numKsamples,reshape(data(x,y,:),1,numKsamples))
                end
            end
        end
    
    %y phase encoding direction
    case 'y'
        for x = 1:xdim;
            for y = 1:ydim;
                if (ismember(x,brain) & ~ismember(y,brain));
                    plot(1:numKsamples,reshape(data(x,y,:),1,numKsamples))
                end
            end
        end
 end


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD STUFF, keeping just in case %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% doCorrectedPhaseRecon %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% function realignedKsp = doCorrectedPhaseRecon(xdim,ydim,numKsamples,originalImage,sampleTimeMs,encodingDirection);
% 
% unalignedKsp = zeros(xdim,xdim,numKsamples);
% realignedKsp = zeros(xdim,xdim,numKsamples);
% 
% % get the kspace image the way they did originally
% for row = 1:xdim
%     %shift the signals a random amount
%     for sample = 1:numKsamples
%         imageKsp = fftshift(fft2(originalImage{row}(:,:,sampleTimeMs*(1+(sample-1)))));      
%         switch encodingDirection
%             case 'x'
%                 kspLine = imageKsp(row,:);
%                 unalignedKsp(row,:,sample) = kspLine;
%                 
%             case 'y'
%                 kspLine = imageKsp(:,row);
%                 unalignedKsp(:,row,sample) = kspLine;
%         end
%     end
%     realignedLine = phaseAlignKsp(unalignedKsp,row,encodingDirection,numKsamples);
%     realignedKsp(row,:,:) = realignedLine;
% end


%%%%%%%%%%%%%%%%%%%
%% phaseAlignKsp %%
%%%%%%%%%%%%%%%%%%%
% function realignedLine = phaseAlignKsp(unalignedKsp,row,encodingDirection,numKsamples)
% 
% %get just the single line
% singleLineFxFyT = zeros(size(unalignedKsp));
% singleLineFxFyT(row,:,:) = unalignedKsp(row,:,:);
% 
% 
% %recon an XYT image from the FxFyT matrix
% for sample = 1:numKsamples,
%     singleLineXYT(:,:,sample) = ifft2(singleLineFxFyT(:,:,sample),'symmetric');
% end
% 
% %take the fourier transform in the time dimension
% singleLineXYFt = fft(singleLineXYT,[],3);
% 
% %take the absolute value to throw out temporal phase information
% absSingleLineXYFt = abs(singleLineXYFt);
% %absSingleLineXYFt = singleLineXYFt;
% 
% %go back to XYT with temporal phase thrown out
% nophaseSingleLineXYT = ifft(absSingleLineXYFt,[],3);
% 
% %go back to FxFyT now that you've removed temporal phase
% for sample = 1:numKsamples,
%     nophaseSingleLineFxFyT(:,:,sample) = fft2(nophaseSingleLineXYT(:,:,sample));
% end
% 
% realignedLine = nophaseSingleLineFxFyT(row,:,:);




%%%%%%%%%%%%%%%%%%%%%%%
%% makeOriginalImageHeartbeat %%
%%%%%%%%%%%%%%%%%%%%%%%
% function originalImage = makeOriginalImageHeartbeat(xdim,ydim,noiseLevel,ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz);
% 
% %make the image with a "brain" which will have signals in it
% for lineSample = 1:xdim
%     originalImage{lineSample} = normrnd(0,noiseLevel,xdim,ydim,ntimePointsMs);
% end
% 
% %put a signal of determined frequency but random phase in the "brain" voxels
% brain = (round(xdim/2)-brainSize/2):(round(xdim/2)+brainSize/2); %brain is square here, and xdim=ydim
% 
% for voxelx = 1:xdim;
%     for voxely = 1:ydim;
%         if ismember(voxelx,brain) & ismember(voxely,brain)
%             voxelOffset = round(rand*ntimePointsMs); %makes the voxels different phase from eachother
%             %voxelOffset = 0;
%             parfor lineSample = 1:xdim
%                 originalImage{lineSample}(voxelx,voxely,:) = brainBaseContrast + normpdf(1:ntimePointsMs,300+round(rand*50),40)*8*(rand/2+1);
%             end
% 
%         end
%     end
% end
