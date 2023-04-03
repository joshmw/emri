%load image, make 101x101, take the FT
    im = imread('peppers.png'); im = rgb2gray(im);
    im = im(100:200,200:300);
    im = im2double(im);
    ksp = fftshift(fft2(im));

%show example of conjugate pairs - i think you get complex values if you
%take them on just one conjugate pair, so you need to take the absolute
%value and it looks normal. I accidentally made kspace 1 line smaller and
%for some reason that didn't give me complex values at one point... worth
%asking about?

    figure
    testksp = zeros(101);
    testksp(44,44) = ksp(44,44); testksp(58,58) = ksp(58,58);
    subplot(3,2,1), imshow(abs(testksp)), title('k space')
    subplot(3,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar
    
    testksp = zeros(101);
    testksp(48,53) = ksp(48,53); testksp(54,49) = ksp(54,49);
    subplot(3,2,3), imshow(abs(testksp)), title('k space')
    subplot(3,2,4), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    subplot(3,2,5), imshow(abs(testksp)), title('k space')
    subplot(3,2,6), imagesc(abs(ifft2(testksp))), title('abs value of above image'), colorbar

%changing the numbers changes contrast, luminance and phase - but not
%orientation of frequency. it also returns complex values... maybe this is
%our problem.
    figure
    subplot(3,2,1), imagesc(abs(testksp)), colorbar, title('k space')
    subplot(3,2,2), imagesc(abs(ifft2(testksp))), colorbar, title('image')
    
    %multiply a single point by 2 for example
    testksp(48,53) = ksp(48,53); testksp(54,49) = ksp(54,49)*5;
    subplot(3,2,3), imagesc(abs(testksp)), colorbar, title('k space')
    subplot(3,2,4), imagesc(abs(ifft2(testksp))), colorbar, title('image')

    %multiply the other point by 2
    testksp(48,53) = ksp(48,53)*5; testksp(54,49) = ksp(54,49);
    subplot(3,2,5), imagesc(abs(testksp)), colorbar, title('k space')
    subplot(3,2,6), imagesc(abs(ifft2(testksp))), colorbar, title('image')


%take a line and its conjugate pair instead of just a point
   
    figure
    testksp = zeros(101);
    testksp(49,:) = ksp(49,:); testksp(53,:) = ksp(53,:);
    subplot(3,2,1), imshow(abs(testksp)), title('k space')
    subplot(3,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(101);
    testksp(49,:) = ksp(49,:); %just doing a single line
    subplot(3,2,3), imshow(abs(testksp)), title('k space')
    subplot(3,2,4), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(101);
    testksp(:,49) = ksp(:,49); %do it horizontal, for proof
    subplot(3,2,5), imshow(abs(testksp)), title('k space')
    subplot(3,2,6), imagesc(abs(ifft2(testksp))), title('image'), colorbar



    figure
    testksp = zeros(101);
    testksp(49,:) = ksp(49,:); testksp(53,:) = randperm(ksp(53,:),length(ksp(53,:)));
    subplot(3,2,1), imshow(abs(testksp)), title('k space')
    subplot(3,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar








keyboard
