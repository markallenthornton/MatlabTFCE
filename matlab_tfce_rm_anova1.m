function pcorr = matlab_tfce_rm_anova1(imgs,nperm,H,E,C,dh)
%MATLAB_TFCE_RM_ANOVA1 performs maximal statistic permutation testing
%   pcorr = matlab_tfce_rm_anova1(imgs,nperm,H,E,C,dh) corrects
%   a one-way omnibus ANOVA (means ~=) for multiple comparisons via a 
%   permutation procedure (random condition switching). Maximal f-stats of
%   permuted data are compared with f-stats in the unshuffled data.
%
%   Arguments:
%   imgs -- an m x 1 cell array of 4D (x,y,z,subject) matrices of images
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values.
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation
%
%   Output:
%   pcorr -- wholebrain map of corrected p-values

% calculate matrix size
bsize = size(imgs{1});
nsub = bsize(4);
bsize = bsize(1:3);
levels = size(imgs,1);

% set tranform function
transform = @matlab_tfce_transform;

% calculate implicit mask
cellmeans = cellfun(@(x) sum(x,4),imgs,'UniformOutput',false);
implicitmask = sum(cat(4,cellmeans{:}),4);
implicitmask = ~(implicitmask == 0 | isnan(implicitmask));
nvox = sum(implicitmask(:));

% extract occupied voxels for permutation test
occimgs = NaN(nsub,levels,nvox);
for i = 1:levels
    for s = 1:nsub
        curimg = imgs{i}(:,:,:,s);
        occimgs(s,i,:) = curimg(implicitmask);
    end
end
clear imgs curimg

% calculate means
smeans = mean(occimgs,2);
cmeans = mean(occimgs,1);
gmeans = mean(mean(occimgs));
gmeans_c = repmat(gmeans,[1,levels,1]);
gmeans_s = repmat(gmeans,[nsub,1,1]);

% calculate SS
SSfac = sum((cmeans-gmeans_c).^2,2)*nsub;
SSsub = sum((smeans-gmeans_s).^2,1)*levels;
SStot = NaN([1,1,nvox]);
for v = 1:nvox
    diffs = occimgs(:,:,v)-gmeans(v);
    SStot(v) = sum(diffs(:).^2);
end
SSw = SStot-SSsub;
SSerr = SSw-SSfac;

% calculate MS and F
MSfac = SSfac./(levels-1);
MSerr = SSerr./((levels-1)*(nsub-1));
fvals = MSfac./MSerr;

% perform TFCE
truestat = zeros(bsize);
truestat(implicitmask) = fvals;
tfcestat = transform(truestat,H,E,C,dh);
tfcestat = tfcestat(implicitmask);

% initialize progress indicator
parfor_progress(nperm);
global parworkers

% cycle through permutations
exceedances = zeros(nvox,1);
parfor(p = 1:nperm,parworkers)
    
    % permute labels
    roccimgs = NaN(size(occimgs));
    for s = 1:nsub
        relabeling = randperm(levels)';
        roccimgs(s,relabeling,:) = occimgs(s,:,:);
    end
    
    % calculate permutation statistic
    rcmeans = mean(roccimgs,1);
    rSSfac = sum((rcmeans-gmeans_c).^2,2)*nsub;
    rSSerr = SSw-rSSfac;
    rMSfac = rSSfac./(levels-1);
    rMSerr = rSSerr./((levels-1)*(nsub-1));
    rfvals = rMSfac./rMSerr;
    rbrain = zeros(bsize);
    rbrain(implicitmask) = rfvals;
    rbrain = transform(rbrain,H,E,C,dh);
    rstats = rbrain(implicitmask);
    
    % compare maxima to true statistic and increment as appropriate
    curexceeds = max(rstats) >= tfcestat;
    exceedances = exceedances + curexceeds;
    
    % update progress indicator (only does so 1 in 5 to minimize overhead)
    if ~randi([0 4]);
        parfor_progress;
    end
    
end

% create corrected p-value image
corrected = exceedances./nperm;
pcorr = ones(bsize);
pcorr(implicitmask) = corrected;

end

