% do the same as above, but shifting the phase of the voxels
function phaseUnlockedKsp = inImageReconPhaseUnlocked(params, inImage);

%can sample 2 lines at once if you'd like to make kspace conjugate symmetric.
%this ensures a REAL IMAGE. effectively the same as half-fourier sampling.
if params.buildConjugateLines

    for row = 1:(floor(params.xdim/2)+1)
        for sample = 1:params.numKsamples
            %takes the temporal order of the images
            imageKsp = fftshift(fft2(inImage{row}(:,:,1+params.sampleTimeMs*(sample-1))));      
            switch params.encodingDirection
                case 'x'
                    kspLine = imageKsp(row,:);
                    phaseUnlockedKsp(row,:,sample) = kspLine;
                    kspLine2 = imageKsp(params.xdim-row+1,:);
                    phaseUnlockedKsp(params.xdim-row+1,:,sample) = kspLine2;
                case 'y'
                    kspLine = imageKsp(:,row);
                    phaseUnlockedKsp(:,row,sample) = kspLine;
                    kspLine2 = imageKsp(:,params.xdim-row+1);
                    phaseUnlockedKsp(:,params.xdim-row+1,sample) = kspLine2;
            end
        end
    end
    
else
    %sampling one line at a time, no conjugate pairing.
    for row = 1:(params.xdim)
        for sample = 1:params.numKsamples
            %takes the temporal order of the images
            imageKsp = fftshift(fft2(inImage{row}(:,:,1+params.sampleTimeMs*(sample-1))));      
            switch params.encodingDirection
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