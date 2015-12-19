% tfce_artifical_test.m

%% test 1: tfce on single cubic cluster
img = zeros([80 80 80]);
img(21:60,21:60,21:60) = randn([40 40 40]);
img(31:50,31:50,31:50) = img(31:50,31:50,31:50) + 1;
img = smooth3(img);
tic
tfced = tfce_transform(img,2,.5,6,10);
toc

figure
subplot(2,2,1)
imagesc(img(:,:,40));
subplot(2,2,2)
plot(img(:,40,40));

subplot(2,2,3)
imagesc(tfced(:,:,40));
subplot(2,2,4)
plot(tfced(:,40,40));

%% test 2: tfce transform on 1-D curve
curve = 100-(-10:1:10).^2;
img = zeros(100,1);
img(10:30) = curve;
img(70:90) = curve*2;
img = img + randn(100,1)*30;
tfced = tfce_transform(img);

figure
subplot(2,1,1)
plot(img)
subplot(2,1,2)
plot(tfced)

%% test 3: permutation on small tfced brain
% simulate images
img = zeros([40 40 20]);
img(16:25,16:25,9:12) = 2;
imgs = repmat(img,[1 1 1 15])+randn([40,40,20 15]);

% perform permutation
tic;pcorr= tfce_permutation(imgs,100); toc

% visualize slice
figure
subplot(1,2,1)
imagesc(pcorr(:,:,10),[0 1]);
subplot(1,2,2)
imagesc(pcorr(:,:,10)<.05,[0 1]);
%% compare with & without stepdown
tic;pcorr_tfced = tfce_permutation(tfced,1000); toc
tic;pcorr_tfced_nostepdown = tfce_permutation_nostepdown(tfced,1000); toc

figure
subplot(2,2,1)
imagesc(pcorr_tfced_nostepdown(:,:,10),[0 1]);
subplot(2,2,2)
imagesc(pcorr_tfced_nostepdown(:,:,10)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr_tfced(:,:,10),[0 1]);
subplot(2,2,4)
imagesc(pcorr_tfced(:,:,10)<.05,[0 1]);

figure
plot(pcorr_tfced_nostepdown(:),pcorr_tfced(:),'.')
hold on
plot(0:1,0:1,'r')

%% compare with & without stepdown (uneven effect sizes)
img = zeros([40 40 20]);
img(21:25,21:25,9:12) = 2;
img(16:20,16:20,9:12) = 5;
imgs = repmat(img,[1 1 1 15])+randn([40,40,20 15]);

tic;pcorr = tfce_permutation(imgs,1000); toc
tic;pcorr_nostepdown = tfce_permutation_nostepdown(imgs,1000); toc

figure
subplot(2,2,1)
imagesc(pcorr_nostepdown(:,:,10),[0 1]);
subplot(2,2,2)
imagesc(pcorr_nostepdown(:,:,10)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr(:,:,10),[0 1]);
subplot(2,2,4)
imagesc(pcorr(:,:,10)<.05,[0 1]);

figure
imagesc(pcorr(:,:,10)-pcorr_nostepdown(:,:,10));colorbar

figure
plot(pcorr_nostepdown(:),pcorr(:),'.')
hold on
plot(0:1,0:1,'r')

%% compare one- and two-tailed tests
img = zeros([40 40 20]);
img(21:25,21:25,9:12) = 2;
img(16:20,16:20,9:12) = 5;
imgs = repmat(img,[1 1 1 15])+randn([40,40,20 15]);

tic;pcorr1 = tfce_permutation(imgs,1000); toc
tic;pcorr2 = tfce_permutation_twotailed(imgs,1000); toc

figure
subplot(2,2,1)
imagesc(pcorr1(:,:,10),[0 1]);
subplot(2,2,2)
imagesc(pcorr1(:,:,10)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr2(:,:,10),[0 1]);
subplot(2,2,4)
imagesc(pcorr2(:,:,10)<.05,[0 1]);

figure
plot(pcorr1(:),pcorr2(:),'.')
hold on
plot(0:1,0:1,'r')

%% two-tailed activation and deactivation
img = zeros([40 40 20]);
img(21:25,21:25,9:12) = -2;
img(16:20,16:20,9:12) = 2;
imgs = repmat(img,[1 1 1 15])+randn([40,40,20 15]);
tic

tic;[pcorr_pos,pcorr_neg] = tfce_permutation_twotailed(imgs,1000); toc

figure
subplot(2,2,1)
imagesc(pcorr_pos(:,:,10),[0 1]);
subplot(2,2,2)
imagesc(pcorr_pos(:,:,10)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr_neg(:,:,10),[0 1]);
subplot(2,2,4)
imagesc(pcorr_neg(:,:,10)<.05,[0 1]);

