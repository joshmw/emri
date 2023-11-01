

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD STUFF, keeping just in case %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% doCorrectedPhaseRecon %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% function realignedKsp = doCorrectedPhaseRecon(xdim,ydim,numKsamples,originalImage,params.sampleTimeMs,encodingDirection);
% 
% unalignedKsp = zeros(xdim,xdim,numKsamples);
% realignedKsp = zeros(xdim,xdim,numKsamples);
% 
% % get the kspace image the way they did originally
% for row = 1:xdim
%     %shift the signals a random amount
%     for sample = 1:numKsamples
%         imageKsp = fftshift(fft2(originalImage{row}(:,:,params.sampleTimeMs*(1+(sample-1)))));      
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
% function originalImage = makeOriginalImageHeartbeat(xdim,ydim,noiseLevel,params.ntimePointsMs,brainSize,brainBaseContrast,amplitude,hz);
% 
% %make the image with a "brain" which will have signals in it
% for lineSample = 1:xdim
%     originalImage{lineSample} = normrnd(0,noiseLevel,xdim,ydim,params.ntimePointsMs);
% end
% 
% %put a signal of determined frequency but random phase in the "brain" voxels
% brain = (round(xdim/2)-brainSize/2):(round(xdim/2)+brainSize/2); %brain is square here, and xdim=ydim
% 
% for voxelx = 1:xdim;
%     for voxely = 1:ydim;
%         if ismember(voxelx,brain) & ismember(voxely,brain)
%             voxelOffset = round(rand*params.ntimePointsMs); %makes the voxels different phase from eachother
%             %voxelOffset = 0;
%             parfor lineSample = 1:xdim
%                 originalImage{lineSample}(voxelx,voxely,:) = brainBaseContrast + normpdf(1:params.ntimePointsMs,300+round(rand*50),40)*8*(rand/2+1);
%             end
% 
%         end
%     end
% end
