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
img(16:25,16:25,9:12) = 1;
imgs = repmat(img,[1 1 1 15])+randn([40,40,20 15]);

% perform tfce
tic
tfced = NaN(size(imgs));
for s = 1:15
    tfced(:,:,:,s) = tfce_transform(imgs(:,:,:,s),2,.5,6,100);
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
%% permutation histograms
x = randn(30,3);
x(:,1) = x(:,1)+5;
x(:,2) = x(:,2)+1;
mvec = zeros(1000,3);
for i = 1:1000
    relabeling = randsample([-1 1],30,'true');
    rx = x;
    for j = 1:30
        if relabeling(j) == -1
          rx(j,:) = -rx(j,:);
        end
    end
    %mvec(i,:) = mean(rx);
    mvec(i,:) = mean(rx)./(std(rx)/sqrt(30));
end
subplot(3,1,1)
hist(mvec(:,1))
xlim([-3 3])
subplot(3,1,2)
hist(mvec(:,2))
xlim([-3 3])
subplot(3,1,3)
hist(mvec(:,3))
xlim([-3 3])

mean(mvec(:)>mean(x(:,2)))
protected = mvec(:,2:3);
mean(protected(:)>mean(x(:,2)))
%% testing max timing
vals = randn(1000,1);
tic
for i = 1:999
    max(vals(1:2));
end
toc
tic
max(vals);
toc
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
tic;[pcorr_pos pcor_neg]= tfce_permutation_twotailed(imgs,1000); toc

figure
subplot(2,2,1)
imagesc(pcorr1(:,:,10),[0 1]);
subplot(2,2,2)
imagesc(pcorr1(:,:,10)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr_pos(:,:,10),[0 1]);
subplot(2,2,4)
imagesc(pcorr_pos(:,:,10)<.05,[0 1]);

figure
plot(pcorr1(:),pcorr_pos(:),'.')
hold on
plot(0:1,0:1,'r')

%% two-tailed activation and deactivation
img = zeros([40 40 20]);
img(21:25,21:25,9:12) = -1;
img(16:20,16:20,9:12) = 1;
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

%% correlation test case 1
covariate = zeros([1,1,1,30]);
covariate(:) = (1:30)-30;
imgs = randn([20 20 10 30])*10;
imgs(6:10,6:10,3:5,1:30) = repmat(covariate,[5,5,3,1]) + imgs(6:10,6:10,3:5,1:30);

% perform tfce (must be two-tailed even if test one-tailed)
tic
tfced = NaN(size(imgs));
for s = 1:30
    tfced(:,:,:,s) = tfce_transform_twotailed(imgs(:,:,:,s),2,.5,6,100);
end
toc

tic;pcorr_imgs = tfce_correlation(imgs,covariate,1000);toc
tic;pcorr_tfced = tfce_correlation(tfced,covariate,1000);toc

figure
subplot(2,2,1)
imagesc(pcorr_imgs(:,:,4),[0 1]);
subplot(2,2,2)
imagesc(pcorr_imgs(:,:,4)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr_tfced(:,:,4),[0 1]);
subplot(2,2,4)
imagesc(pcorr_tfced(:,:,4)<.05,[0 1]);

figure
plot(pcorr_imgs(:),pcorr_tfced(:),'.')
hold on
plot(0:1,0:1,'r')

%% correlation test case 2 (two-tailed)
covariate = zeros([1,1,1,30]);
covariate(:) = (1:30)-30;
imgs = randn([20 20 10 30])*10;
imgs(6:10,6:10,3:5,1:30) = repmat(covariate,[5,5,3,1]) + imgs(6:10,6:10,3:5,1:30);
imgs(11:15,11:15,3:5,1:30) = repmat(-covariate,[5,5,3,1]) + imgs(11:15,11:15,3:5,1:30);

% perform tfce
tic
tfced = NaN(size(imgs));
for s = 1:30
    tfced(:,:,:,s) = tfce_transform_twotailed(imgs(:,:,:,s),2,.5,6,100);
end
toc

tic;[pcorr_pos,pcorr_neg] = tfce_correlation_twotailed(tfced,covariate,1000);toc

figure
subplot(2,2,1)
imagesc(pcorr_pos(:,:,4),[0 1]);
subplot(2,2,2)
imagesc(pcorr_pos(:,:,4)<.05,[0 1]);
subplot(2,2,3)
imagesc(pcorr_neg(:,:,4),[0 1]);
subplot(2,2,4)
imagesc(pcorr_neg(:,:,4)<.05,[0 1]);

%% independent samples test
img = zeros([40 40 20]);
img(21:25,21:25,9:12) = 2;
imgs1 = repmat(img,[1 1 1 15])+randn([40,40,20 15]);
img = zeros([40 40 20]);
img(16:20,16:20,9:12) = 2;
imgs2 = repmat(img,[1 1 1 20])+randn([40,40,20 20]);

% perform tfce
tic
tfced1 = NaN(size(imgs1));
for i = 1:15
    tfced1(:,:,:,i) = tfce_transform_twotailed(imgs1(:,:,:,i),2,.5,6,100);
end
tfced2 = NaN(size(imgs2));
for i = 1:20
    tfced2(:,:,:,i) = tfce_transform_twotailed(imgs2(:,:,:,i),2,.5,6,100);
end
toc

% perform permutation tests
tic;pcorr = tfce_permutation_independent(tfced1,tfced2,100);toc
tic;[pcorr_pos,pcorr_neg] = tfce_permutation_independent_twotailed(tfced1,tfced2,100);toc

% visualize results
subplot(3,2,1)
imagesc(pcorr(:,:,10),[0,1])
subplot(3,2,2)
imagesc(pcorr(:,:,10)<.05,[0,1])
subplot(3,2,3)
imagesc(pcorr_pos(:,:,10),[0,1])
subplot(3,2,4)
imagesc(pcorr_pos(:,:,10)<.05,[0,1])
subplot(3,2,5)
imagesc(pcorr_neg(:,:,10),[0,1])
subplot(3,2,6)
imagesc(pcorr_neg(:,:,10)<.05,[0,1])
