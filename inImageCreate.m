% function inImage = inImageCreate(imType,xdim,ydim,noiseLevel,ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz,offsetVoxels,phaseShiftBetweenLines)
function inImage = inImageCreate(imType,eP)
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

%% Parse input parameters
%{
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
%}

inImage = cell(eP.xdim,1);

% Make the image noise samples
for lineSample = 1:eP.xdim
    % Each cell is a noisy video.
    inImage{lineSample} = normrnd(0,eP.noiseLevel,eP.xdim,eP.ydim,eP.ntimePointsMs);
end

switch imType
    case {'sine','harmonic'}

        % Shifts the time series for a specific x-line (each x-line is
        % a cell)
        if phaseShiftBetweenLines
            lineOffsetMs = round(rand(1,xdim)*ntimePointsMs);
        end

        % Determine the brain voxel locations
        % Brain is square for debugging
        brain = (round(xdim/2)-brainSize/2):(round(xdim/2)+brainSize/2); 

        for voxelx = 1:xdim
            for voxely = 1:ydim
                if ismember(voxelx,brain) && ismember(voxely,brain)

                    % Shifts the time series across all cells for a
                    % specific voxel
                    if offsetVoxels
                        % makes the voxels different phase from each other
                        voxelOffset = round(rand(1,1) * ntimePointsMs); 
                    else
                        voxelOffset = 0;
                    end

                    % make the lines
                    if ~phaseShiftBetweenLines
                        % The time points depend only on the line you
                        % are sampling, without any additional offset.                        
                        parfor lineSample = 1:xdim
                            t = ((1+(lineSample-1)*ntimePointsMs):(lineSample*ntimePointsMs));
                            inImage{lineSample}(voxelx,voxely,:) = ...
                                brainBaseContrast + amplitude*cos(t*2*pi/1000*hz+voxelOffset);
                        end
                    else
                        parfor lineSample = 1:xdim
                            % The time points depend on a random
                            % offset the is determined by the line you
                            % are sampling.
                            t = ((1+lineOffsetMs(lineSample)):(lineOffsetMs(lineSample))+ntimePointsMs);
                            inImage{lineSample}(voxelx,voxely,:) = ...
                                brainBaseContrast + amplitude*cos(t*2*pi/1000*hz + voxelOffset);
                        end
                    end
                end
            end
        end
    otherwise
        error('Unknown inImage type %s\n',imType);
end

end

