%load image, make 101x101, take the FT
    im = imread('peppers.png'); im = rgb2gray(im);
    im = im(100:300,200:400);
    im = im2double(im); %im = im-mean(im(:));
    ksp = fftshift(fft2(im));

    figure
    subplot(1,2,1), imshow(im); title('original Image')
    subplot(1,2,2), imshow(abs(ksp)); title('abs(fftshift(fft2(image)))')
    sgtitle('Original image and abs of fftshift of fft2')

%show example of conjugate pairs - i think you get complex values if you
%take them on just one conjugate pair, so you need to take the absolute
%value and it looks normal. I accidentally made kspace 1 line smaller and
%for some reason that didn't give me complex values at one point... worth
%asking about?

keyboard

    figure
    testksp = zeros(201);
    testksp(94,94) = ksp(94,94); testksp(108,108) = ksp(108,108);
    subplot(2,2,1), imshow(abs(testksp)), title('k space')
    subplot(2,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar
    
    testksp = zeros(201);
    testksp(98,103) = ksp(98,103); testksp(104,99) = ksp(104,99);
    subplot(2,2,3), imshow(abs(testksp)), title('k space')
    subplot(2,2,4), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    %subplot(3,2,5), imshow(abs(testksp)), title('k space')
    %subplot(3,2,6), imagesc(abs(ifft2(testksp))), title('abs value of above image'), colorbar

    sgtitle('ifft2 of single points + conjugate pairs')

%changing the numbers changes contrast, luminance and phase - but not
%orientation of frequency. it also returns complex values... maybe this is
%our problem.
keyboard

    figure
    subplot(3,2,1), imagesc(abs(testksp)), colorbar, title('k space (top right point x5)')
    subplot(3,2,2), imagesc(abs(ifft2(testksp))), colorbar, title('image')
    
    %multiply a single point by 2 for example
    testksp(98,103) = ksp(98,103); testksp(104,99) = ksp(104,99)*5;
    subplot(3,2,3), imagesc(abs(testksp)), colorbar, title('k space (bottom left point x5)')
    subplot(3,2,4), imagesc(abs(ifft2(testksp))), colorbar, title('image')

    %multiply the other point by 2
    testksp(98,103) = ksp(98,103); testksp(104,99) = ksp(round(rand*200),round(rand*200));
    subplot(3,2,5), imagesc(abs(testksp)), colorbar, title('k space - cojugate point is random from elsewhere in image')
    subplot(3,2,6), imagesc(abs(ifft2(testksp))), colorbar, title('image')

    sgtitle('ifft2 of single points + scaled/swapped conjugate pairs')


%take a line and its conjugate pair instead of just a point
keyboard   

    figure
    testksp = zeros(201);
    testksp(99,:) = ksp(99,:); testksp(103,:) = ksp(103,:);
    subplot(3,2,1), imshow(abs(testksp)), title('k space, paired line')
    subplot(3,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(201);
    testksp(99,:) = ksp(99,:); %just doing a single line
    subplot(3,2,3), imshow(abs(testksp)), title('k space, single line')
    subplot(3,2,4), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(201);
    testksp(:,99) = ksp(:,99); %do it horizontal, for proof
    subplot(3,2,5), imshow(abs(testksp)), title('k space')
    subplot(3,2,6), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    sgtitle('ifft2 of lines, with/without pairs')

keyboard
    figure
    testksp = zeros(201);
    testksp(99,:) = ksp(99,:); testksp(103,1:10) = ksp(103,1:10);
    subplot(4,2,1), imshow(abs(testksp)), title('k space, section of conjugate line')
    subplot(4,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(201);
    testksp(99,:) = ksp(99,:); testksp(103,96:106) = ksp(103,96:106);
    subplot(4,2,3), imshow(abs(testksp)), title('k space, section of conjugate line')
    subplot(4,2,4), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(201);
    testksp(99,:) = ksp(99,:); testksp(103,:) = ksp(140,:);
    subplot(4,2,5), imshow(abs(testksp)), title('k space, conjugate line is different random line')
    subplot(4,2,6), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    testksp = zeros(201);
    testksp(:,99) = ksp(:,99); testksp(:,103) = ksp(:,130);
    subplot(4,2,7), imshow(abs(testksp)), title('k space, conjugate line is different random line')
    subplot(4,2,8), imagesc(abs(ifft2(testksp))), title('image'), colorbar
   
   sgtitle('ifft2 of different combinations of lines')
















   %%%% just testing things here %%%%%


keyboard

testksp = zeros(201);
for i = 1:201;
    noisyIm = im;
    noisyIm(i,51:151) = noisyIm(i,51:151)+normrnd(0,.1,1,101);
    noisyKsp = fftshift(fft2(noisyIm));
    testksp(i,:) = noisyKsp(i,:);
end
figure,imshow(abs(ifft2(testksp)))



figure,imshow(abs(ifft2(testksp)-im))














figure
randKsp = ksp(:,randperm(size(ksp, 1)))
subplot(1,2,1), imshow(abs(randKsp)), title('k space, conjugate line is different random line')
subplot(1,2,2), imagesc(abs(ifft2(randKsp))), title('image'), colorbar






%checking if adding ifft2 of non-conjugate paired lines is linear (it is)
    figure
    testksp = zeros(201);
    testksp(99,:) = ksp(99,:); testksp(103,96:106) = ksp(103,96:106);
    subplot(2,2,1), imshow(abs(testksp)), title('k space, section of conjugate line')
    subplot(2,2,2), imagesc(abs(ifft2(testksp))), title('image'), colorbar

    
    testksp1 = zeros(201); testksp2 = zeros(201);
    testksp1(99,:) = ksp(99,:);
    testksp2(103,96:106) = ksp(103,96:106)
    testImage = ifft2(testksp1) + ifft2(testksp2);
    subplot(2,2,4), imagesc(abs(testImage))









''