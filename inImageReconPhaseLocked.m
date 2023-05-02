function ksp = inImageReconPhaseLocked(params,inImage)
% Purpose:
%  Samples an inImage the way described in Toi et al. Returns a frequency movie analogous to k-space.
%  Assumes *perfect stationarity* of the underlying image between lines by only sampling the same movie (first cell of inImage) for each line.
%
% Usage: 
%  params = emriParams;
%  inImage = inImageCreate(params,signalType)
%  ksp = inImageReconPhaseLocked(params,inImage)
%
% Output:
%  ksp: xdim x xdim x ntimePointsMS/sampleTimeMs
%
% Required:
%   inImage: cell array of XYT movies from inImageCreate. 
%   params: Set during emriSimul.
%    'xdim': size of the image (x by x)
%    'ntimePointsMs': length of one line acquisition
%    'sampleTimeMs': TR length. Takes a sample (1 line) every n miliseconds. i.e. nTimePointsMs = 500 and sampleTimeMS = 5 -> 100 samples of each line.
%    'numKsamples': Number of measurements. Time series length/TR length (ntimePointsMS/sampleTimeMs).



%initialize kspace for speed
ksp = zeros(params.xdim,params.ydim,params.numKsamples);

% sample each line (row or column) at each sample timepoint (TR length)
for row = 1:params.xdim
    %for each TR, sample the line in kspace and add to the kspace image
    for sample = 1:params.numKsamples  
        imageKsp = fftshift(fft2(inImage{1}(:,:,1+params.sampleTimeMs*(sample-1))));      
        switch params.encodingDirection
            case 'x'
                kspLine = imageKsp(row,:);
                ksp(row,:,sample) = kspLine;
            case 'y'
                kspLine = imageKsp(:,row);
                ksp(:,row,sample) = kspLine;
        end
    end
end