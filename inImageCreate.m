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
    inImage{lineSample} = normrnd(0,params.noiseLevel,params.xdim,params.ydim,params.ntimePointsMs);
end

switch params.signalType
    case {'sine','harmonic'}

        % Shifts the time series for a specific x-line (each x-line is
        % a cell)
        if params.phaseShiftBetweenLines
            lineOffsetMs = round(rand(1,params.xdim)*params.ntimePointsMs);
        end

        % Determine the brain voxel locations
        % Brain is square for debugging
        brain = (round(params.xdim/2)-params.brainSize/2):(round(params.xdim/2)+params.brainSize/2); 

        for voxelx = 1:params.xdim
            for voxely = 1:params.ydim
                if ismember(voxelx,brain) && ismember(voxely,brain)

                    % Shifts the time series across all cells for a specific voxel
                    if params.offsetIndividualVoxels
                        % makes the voxels different phase from each other
                        voxelOffset = round(rand(1,1) * params.ntimePointsMs); 
                    else
                        voxelOffset = 0;
                    end

                    % make the time series for each cell
                    if ~params.phaseShiftBetweenLines
                        % The time points depend only on the line you are sampling, without any additional offset.                        
                        parfor lineSample = 1:params.xdim
                            t = ((1+(lineSample-1)*params.ntimePointsMs):(lineSample*params.ntimePointsMs));
                            inImage{lineSample}(voxelx,voxely,:) = ...
                                params.brainBaseContrast + params.amplitude*cos(t*2*pi/1000*params.hz+voxelOffset);
                        end
                    else
                        % The time points depend on a random offset that is determined by the line you are sampling.
                        parfor lineSample = 1:params.xdim
                            t = ((1+lineOffsetMs(lineSample)):(lineOffsetMs(lineSample))+params.ntimePointsMs);
                            inImage{lineSample}(voxelx,voxely,:) = ...
                                params.brainBaseContrast + params.amplitude*cos(t*2*pi/1000*params.hz + voxelOffset);
                        end
                    end
                end
            end
        end
    otherwise
        error('Unknown inImage type %s\n',imType);
end

