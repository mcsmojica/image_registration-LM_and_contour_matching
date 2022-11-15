% removes noise from hands data and smoothes the pixelated edges so that LM
% detection can be done on binary images
dataR = double(imread('hands-R.jpg'));
dataT = double(imread('hands-T.jpg'));
dataR = imresize(dataR,4);
dataT = imresize(dataT,4);

% remove noise
R_Kaverage = filter2(fspecial('average',3),dataR)/max(dataR(:));
T_Kaverage = filter2(fspecial('average',3),dataT)/max(dataT(:));

% binarize images
R_bw = double(imbinarize(R_Kaverage));
T_bw = double(imbinarize(T_Kaverage));

% smooth the pixelated edges
windowSize = 10;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(R_bw), kernel, 'same');
R_bw = blurryImage > 0.5; % 
blurryImage = conv2(single(T_bw), kernel, 'same');
T_bw = blurryImage > 0.5; % 
dataR_orig = dataR; dataT_orig = dataT;
dataR = R_bw; dataT = T_bw;

clearvars -except dataR dataT R_bw T_bw dataR_orig dataT_orig