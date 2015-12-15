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
tic
img = zeros([40 40 20]);
img(11:30,11:30,6:15) = randn([20,20,10]);
img(16:25,16:25,9:12) = img(16:25,16:25,9:12)+1;
imgs = repmat(img,[1 1 1 15])+randn([40,40,20 15]);
toc

% perform tfce
tic
tfced = NaN(size(imgs));
for s = 1:15
    tfced(:,:,:,s) = tfce_transform(imgs(:,:,:,s),2,.5,6,10);
end
toc

% perform permutation
tic;pcorr_voxelwise = tfce_permutation(imgs,100); toc
tic;pcorr_tfced = tfce_permutation(tfced,100); toc

% visualize slice
figure
subplot(2,2,1)
imagesc(pcorr_voxelwise(:,:,10),[0 1]);
subplot(2,2,2)
imagesc(pcorr_voxelwise(:,:,10)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr_tfced(:,:,10),[0 1]);
subplot(2,2,4)
imagesc(pcorr_tfced(:,:,10)<.05,[0 1]);