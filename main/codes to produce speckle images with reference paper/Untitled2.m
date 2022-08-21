gamma = 5*5; % M to N ratio. set this square number. only for #visualization. (related to line 29)
SNR = 100; % linear signal-to-noise ratio (SNR) to mimic practical noisy situations.
%% Incident field load
%load('example_inputField.mat'); % load variable X, size(X) = [N, 1] complex double matrix.
X=TestTest1{3,1}(:,:);
N = numel(X); % set this square number. only for #visualization.
M = gamma*N;
X = X ./ norm(X); % normalize X for simplicity
X2d = reshape(X,[sqrt(N),sqrt(N)]); % 2D reshape. only for #visualization.
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
padarray(fftshift(fft2(ifftshift(Y_prac2d))),[2*sqrt(M),2*sqrt(M)]))))).^2; % oversample Iyy. only for #visualization.
%% Field retrieval using SSM method
disp('initial guess start!'); tic;
sigma = mean(abs(TM_prac).^2,1); % normalization factor, sigma (see Eq. (2)).
Zmatrix = (TM_prac'*(repmat(Iyy-mean(Iyy),[1,N]).*TM_prac))./(sigma'*sigma); % Z matrix calculation.
[Xretv0,~] = eigs(Zmatrix,1); % Field retrieval by selecting the eigenvector of the largest eigenvalue.
disp(['initial guess end! ... Elapsed time: ', num2str(toc),'s']);
Xcorr0 = X'*Xretv0/norm(X)/norm(Xretv0); % calculate correlation between X and Xretv0.
Xretv0 = Xretv0*exp(-1i*angle(Xcorr0)); % global phase correction. only for #visualization.
Xretv02d = reshape(Xretv0,[sqrt(N),sqrt(N)]); % 2D reshape. only for #visualization.
%% #Visualization
figure(1),
set(1,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(243),imshow(abs(X2d),[min(abs(X)),max(abs(X))]),axis image; title('incident amplitude'); colorbar
subplot(244),imshow(angle(X2d)),axis image; title('incident wavefront'); colorbar
Xretv02d1 = imresize(abs(Xretv02d), [64 64]);
Xretv02d2 = imresize(Xretv02d1, [64*64 1]);
subplot(247),imshow(Xretv02d1,[min(Xretv02d2),max(Xretv02d2)]),axis on;
%subplot(247),imshow(abs(Xretv02d),[min(abs(X)),max(abs(X))]),axis image; title('retrieved amplitude'); colorbar
subplot(248),imshow(angle(Xretv02d)),axis image; title(['retrieved wavefront, Corr. = ', num2str(abs(Xcorr0))]); colorbar
subplot(2,4,[1,2,5,6]), imshow(Iyy2d_osampd,[0,max(Iyy2d_osampd(:))*0.5]), axis image; title('measured intensity speckle'); 
colormap(subplot(2,4,[1,2,5,6]),'gray'); 
pause(0.1);
%% (optional) Modified Gerchberg-Saxton algorithm
iterMax = 1000; % max iteration number to prevent infinite loop.
iterTor = 10^-5; % convergence criteria.
disp('inverse TM calculation start!'); tic;
[U,S,V]=svd(TM_prac,'econ'); % economic svd for inverse TM calculation.
sv = diag(S);
TMinv_prac = V*(repmat(1./sv,[1,M]).*U'); % inverse TM calculation.
clearvars U S V
disp(['inverse TM calculation end! ... Elapsed time: ', num2str(toc),'s']);
disp('Modified Gerchberg-Saxton iteration start!'); tic;
Xcorr_iter = zeros(iterMax,1); % archives correlation changes during iteration.
Xiter = Xretv0; % set retrieved field for the initial guess.
Xcorr_iter(1) = Xcorr0;
for itr = 1:1:iterMax
Yiter = TM_prac * Xiter; % transform incident field into corresponding diffused field space.
Yiter = sqrt(Iyy).*exp(1i*angle(Yiter)); % replace amplitude part with measured intensity.
Xiter = TMinv_prac*Yiter; % return to the incident field space.
Xcorr_iter(itr+1) = X'*Xiter/norm(Xiter); % correlation calculation.
if abs(Xcorr_iter(itr+1)) - abs(Xcorr_iter(itr)) < iterTor % iteration ends when the criterion is satisfied.
break;
end
end
disp(['Modified Gerchberg-Saxton iteration end! ... Elapsed time: ', num2str(toc),'s']);
Xiter = Xiter*exp(-1i*angle(Xcorr_iter(itr))); % global phase correction. only for #visualization.
Xiter2d = reshape(Xiter,[sqrt(N),sqrt(N)])./ norm(Xiter); % 2D reshape. only for #visualization.
%% #Visualization
figure(2),
set(2,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(241),imshow(abs(X2d),[min(abs(X)),max(abs(X))]),axis image; title('incident amplitude'); colorbar
subplot(242),imshow(angle(X2d)),axis image; title('incident wavefront'); colorbar
Xiter2d1 = imresize(abs(Xiter2d), [64 64]);
Xiter2d2 = imresize(Xiter2d1, [64*64 1]);
subplot(245),imshow(Xiter2d1,[min(Xiter2d2),max(Xiter2d2)]);
%subplot(245),imshow(abs(Xiter2d),[min(abs(X)),max(abs(X))]),axis image; title('retrieved amplitude'); colorbar
subplot(246),imshow(angle(Xiter2d)),axis image; title(['retrieved wavefront, Corr. = ', num2str(abs(Xcorr_iter(itr+1)))]); colorbar
subplot(2,4,[3,4,7,8]), plot(abs(Xcorr_iter(1:itr+1))), title(num2str(itr))
pause(0.1);