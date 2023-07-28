function inImage = inImageCreate(params)
%
% The image is represented as a noisy video for the duration of each
% scan line (xdim). These are stored as a cell array.
%
% Each cell represents the video sequency that could arise as you can
% one of the x-lines. The video in the cell is a matrix of size
%
%   xdim x ydim x ntimePointsMs
%
% See also
%

% Examples:
%{
 eP = emriParams;
 inImage = inImageCreate('harmonic',eP)
%}

% Make the image noise samples
for lineSample = 1:params.xdim
    % Each cell is a noisy video.
    inImage{lineSample} = zeros(params.xdim,params.ydim,params.ntimePointsMs);
    inImageNoise{lineSample} = normrnd(0,params.noiseLevel,params.xdim,params.ydim,params.ntimePointsMs);
end
%definethe brain voxels
brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2); 
%set amp of signal to 0 if not adding
if ~params.addSignal, params.signalAmp = 0; end


%% add the signal.
switch params.signalType
    
    %if you are adding sine modulation)
    case {'sine','harmonic'}

        % Shifts the time series for a specific x-line (each x-line is
        % a cell)
        if params.phaseShiftBetweenLines
            lineOffsetMs = round(rand(1,params.xdim)*params.phaseShiftBetweenLinesMagnitude*params.ntimePointsMs)
        end

        if params.addSignal
            signalOffsetMs = round(rand(1,params.xdim)*params.signalShiftMagnitude)
        else
            signalOffsetMs = zeros(1,params.xdim);
        end
       

        % Determine the brain voxel locations
        % Brain is square for debugging

        for voxelx = 1:params.xdim
            for voxely = 1:params.ydim
                if ismember(voxelx,brain) && ismember(voxely,brain)

                    % Shifts the time series across all cells for a specific voxel
                    if params.offsetIndividualVoxels
                        % makes the voxels different phase from each other
                        voxelOffset = round(rand(1,1) * params.ntimePointsMs * params.offsetIndividualVoxelsMagnitude); 
                    else
                        voxelOffset = 0;
                    end

                    % make the time series for each cell
                    if ~params.phaseShiftBetweenLines
                        % The time points depend only on the line you are sampling, without any additional offset.                        
                        parfor lineSample = 1:params.xdim
                            t = ((1+(lineSample-1)*params.ntimePointsMs):(lineSample*params.ntimePointsMs));
                            inImage{lineSample}(voxelx,voxely,:) = ...        
                                params.brainBaseContrast + params.amplitude*cos(t*2*pi/1000*params.hz+voxelOffset) + ...
                                normpdf((1-signalOffsetMs(lineSample)):(params.ntimePointsMs-signalOffsetMs(lineSample)),params.signalMean,params.signalStd) / ... 
                                max(normpdf((1-signalOffsetMs(lineSample)):(params.ntimePointsMs-signalOffsetMs(lineSample)),params.signalMean,params.signalStd)) * params.signalAmp;
                        end
                    else
                        % The time points depend on a random offset that is determined by the line you are sampling.
                        parfor lineSample = 1:params.xdim
                            t = ((1+lineOffsetMs(lineSample)):(lineOffsetMs(lineSample))+params.ntimePointsMs);
                            inImage{lineSample}(voxelx,voxely,:) = ...
                                params.brainBaseContrast + params.amplitude*cos(t*2*pi/1000*params.hz + voxelOffset) +  ...
                                normpdf((1-signalOffsetMs(lineSample)):(params.ntimePointsMs-signalOffsetMs(lineSample)),params.signalMean,params.signalStd) / ... 
                                max(normpdf((1-signalOffsetMs(lineSample)):(params.ntimePointsMs-signalOffsetMs(lineSample)),params.signalMean,params.signalStd)) * params.signalAmp;
                        end
                    end
                end
            end
        end
    
    %signal other than harmonic 
    otherwise
        error('Unknown inImage type %s\n',imType);
end


    %if you are adding pink noise
    
    if params.addPinkNoise

        if params.addSignal
            signalOffsetMs = round(rand(1,params.xdim)*params.signalShiftMagnitude)
        else
            signalOffsetMs = zeros(1,params.xdim);
        end

        parfor lineSample = 1:params.xdim
            %create pink noise to be used in every voxel for this line
            pinkNoiseLine = squeeze(pinknoise(params.ntimePointsMs))';
            %set the std of the pink noise to be the parameter
            pinkNoiseLine = pinkNoiseLine * params.pinkNoiseStd/std(pinkNoiseLine);
            %loop through the voxels and create time series based on noise
            for voxelx = 1:params.xdim
                for voxely = 1:params.ydim
                    if ismember(voxelx,brain) && ismember(voxely,brain)
                        %add the noise for that line, same across all voxels
                        %add signal if there is one, also same across all voxels
                        inImage{lineSample}(voxelx,voxely,:) = ...
                        inImage{lineSample}(voxelx,voxely,:) + reshape(pinkNoiseLine,[1 1 params.ntimePointsMs]);
                    end
                end
            end
        end    
    end
    



for line = 1:length(inImage);
    inImage{line} = inImage{line}+inImageNoise{line};
end
