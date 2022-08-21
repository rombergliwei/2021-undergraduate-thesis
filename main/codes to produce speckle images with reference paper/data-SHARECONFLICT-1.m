for i=1:2000
gamma = 5*5; % M to N ratio. set this square number. only for #visualization. (related to line 29)
SNR = 100; % linear signal-to-noise ratio (SNR) to mimic practical noisy situations.
%% Incident field load
load('TestTest7.mat'); % load variable X, size(X) = [N, 1] complex double matrix.
X=TestTest7{1,1}(:,:);
N = numel(X); % set this square number. only for #visualization.
M = gamma*N;
X = X ./ norm(X); % normalize X for simplicity
%X2d = reshape(X,[sqrt(N),sqrt(N)]); % 2D reshape. only for #visualization.
%% Transmission matrix (TM) and intensity speckle generation
disp('TM generate start!'); tic;
TM = raylrnd(sqrt(1/2/M),M,N).*exp(1i * random('unif',-pi,+pi-eps('double'),M,N)); % transmission matrix generation.
TM_prac = awgn(TM,SNR,'measured','linear'); % add noise upon preset SNR.
disp(['TM generate end! ... Elapsed time: ', num2str(toc),'s']);
Y = TM * X; % diffused field y generation.
Y_prac = awgn(Y,SNR,'measured','linear'); % add noise upon preset SNR.
Y_prac2d = reshape(Y_prac,[sqrt(M),sqrt(M)]); % 2D reshape. #for visualization.
Iyy = abs(Y_prac).^2; % measured intensity speckle single shot.
Iyy2d_osampd = abs(fftshift(ifft2(ifftshift(...
padarray(fftshift(fft2(ifftshift(Y_prac2d))),[sqrt(M)/10,sqrt(M)/10]))))).^2; % oversample Iyy. only for #visualization.
JIEGUO=mat2gray(Iyy2d_osampd,[0,max(Iyy2d_osampd(:))*0.5]);
B = imresize(JIEGUO, [64 64]);
TestTest8{i,1}(:,:)=B;
%figure;imshow(JIEGUO);
imwrite(B,['test/speckle/',num2str(i),'.jpg']);
end