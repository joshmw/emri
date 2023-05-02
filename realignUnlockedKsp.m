function realignedKsp = kspRealign(params,phaseUnlockedKsp)

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





%%%%%%%%%%%%%%%%%%%%%%%%
%% plotPhaseAlignment %%
%%%%%%%%%%%%%%%%%%%%%%%%
function plotPhaseAlignment(row,params,singleLineFxFyT,singleLineXYT,nophaseSingleLineFxFyT,nophaseSingleLineXYT)
%this only works for 'x' phase encoding direction - could make it work
%otherwise, but I am lazy. will add later.

fig = mlrSmartfig('line image')

%plot a snapshot of the frequency image
subplot(2,4,1)
imagesc(abs(singleLineFxFyT(:,:,1))),colorbar,colormap(gray)
xlabel('Fx'),ylabel('Fy'),title('OG Frequency image (1 timepoint)')

%plot the timecourse of the individual frequencies
subplot(4,4,2), hold on,
for x = 1:params.xdim
    plot(1:params.numKsamples,reshape(real(singleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Real Value Frequency component magnitudes over time (1 row)')

subplot(4,4,6), hold on,
for x = 1:params.xdim
    plot(1:params.numKsamples,reshape(imag(singleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Imaginary Value Frequency component magnitudes over time (1 row)')

%plot a snapshot of the actual image
subplot(2,4,3);
imagesc(singleLineXYT(:,:,1));colorbar,colormap(gray)
xlabel('x'),ylabel('y'),title('image created from single lines')

%plot the timecourse of the image pixels
subplot(2,4,4), hold on
for x = 1:params.xdim
plot(1:params.numKsamples,reshape(singleLineXYT(x,8,:),1,200))
end
xlabel('time'),ylabel('magnitude'),title('pixel values over time (1 column)')

%%then do all of that for the phase-corrected image...
subplot(2,4,5)
imagesc(abs(nophaseSingleLineFxFyT(:,:,1))),colorbar,colormap(gray)
xlabel('Fx'),ylabel('Fy'),title('OG Frequency image (1 timepoint)')

subplot(4,4,10), hold on,
for x = 1:params.xdim
plot(1:params.numKsamples,reshape(real(nophaseSingleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Real Value Frequency component magnitudes over time (1 row)')

subplot(4,4,14), hold on,
for x = 1:params.xdim
plot(1:params.numKsamples,reshape(imag(nophaseSingleLineFxFyT(row,x,:)),1,200))
end
xlabel('time'),ylabel('magnitude of component'),title('Imaginary Value Frequency component magnitudes over time (1 row)')

subplot(2,4,7);
imagesc(nophaseSingleLineXYT(:,:,1));colorbar,colormap(gray)
xlabel('x'),ylabel('y'),title('Realigned image')

subplot(2,4,8), hold on
for x = 1:params.xdim
plot(1:params.numKsamples,reshape(nophaseSingleLineXYT(x,8,:),1,200))
end
xlabel('time'),ylabel('magnitude')
