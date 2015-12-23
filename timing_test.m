%stepdown_tfce('onesample',1,randn([50,50,50,20]),[],[],2);
%%
rng(1);
nsub = 20;
imgs = randn([80,80,80,nsub]);
for i = 1:nsub
    imgs(:,:,:,i) = smooth3(imgs(:,:,:,i),'gaussian',7,1.5);
end
tvals = mean(imgs,4)./(std(imgs,0,4)./sqrt(nsub));
tic;y = stepdown_tfce_transform(tvals,2,.5,6,.1);toc;