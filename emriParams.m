function eParams = emriParams
%  Create a struct with the emri default parameters
%

eParams.xdim = 17; % pixels x
eParams.ydim = 17; % pixels y
eParams.brainSize         = 4;
eParams.encodingDirection = 'x'; %phase encoding direction
eParams.offsetIndividualVoxels = 0; %1 to add a random offset to voxel phase in 'sine' condition
eParams.phaseShiftBetweenLines = 1; % will randomize the phase between lines instead of making continuous 
eParams.noiseLevel = 0;       %std of gaussian noise added to OG image
eParams.ntimePointsMs = 1000; %length of simulation
eParams.hz = 4;                       %signal frequency in "brain" voxels
eParams.amplitude = 1;
eParams.sampleTimeMs = 5; %TR length
eParams.brainBaseContrast = 1;
eParams.buildConjugateLines = 1;

end
