% tfce_artifical_test.m

%% test 1: single square cluster
img = zeros([80 80 80]);
img(21:60,21:60,21:60) = randn([40 40 40]);
img(31:50,31:50,31:50) = img(31:50,31:50,31:50) + 1;
img = smooth3(img);
tfced = tfce_transform(img,2,.5,6,100);

subplot(2,2,1)
imagesc(img(:,:,40));
subplot(2,2,2)
plot(img(:,40,40));

subplot(2,2,3)
imagesc(tfced(:,:,40));
subplot(2,2,4)
plot(tfced(:,40,40));