function realignedKsp = kspRealign(phaseUnlockedKsp,params)

realignedKsp = zeros(size(phaseUnlockedKsp));

for row = 1:(floor(params.xdim/2)+1)

    %take a single pair of FxFyT lines
    singleLineFxFyT = zeros(size(phaseUnlockedKsp));
    switch params.encodingDirection
        case 'x'
            singleLineFxFyT(row,:,:) = phaseUnlockedKsp(row,:,:);
            singleLineFxFyT(params.xdim-row+1,:,:) = phaseUnlockedKsp(params.xdim-row+1,:,:);
        case 'y'
            singleLineFxFyT(:,row,:) = phaseUnlockedKsp(:,row,:);
            singleLineFxFyT(:,params.xdim-row+1,:) = phaseUnlockedKsp(:,params.xdim-row+1,:);
    end

    %recon an XYT image from the single line FxFyT matrix
    for sample = 1:params.numKsamples,
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
    for sample = 1:params.numKsamples,
        nophaseSingleLineFxFyT(:,:,sample) = fftshift(fft2(nophaseSingleLineXYT(:,:,sample)));
    end
    
    for sample = 1:params.numKsamples
        realignedKsp(:,:,sample) = realignedKsp(:,:,sample) + nophaseSingleLineFxFyT(:,:,sample);
        %realignedKsp(row,:,sample) = nophaseSingleLineFxFyT(row,:,sample);
        %realignedKsp(params.xdim-row+1,:,sample) = nophaseSingleLineFxFyT(params.xdim-row+1,:,sample);
    end
  
    if row == floor(params.xdim/2);
        plotPhaseAlignment(row,params,singleLineFxFyT,singleLineXYT,nophaseSingleLineFxFyT,nophaseSingleLineXYT);
    end

end

