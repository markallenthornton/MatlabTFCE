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

% calculate true f-stat image
fvals = NaN(nvox,1);
for v = 1:nvox
    [~,table] = anova_rm(occimgs(:,:,v),'off');
    fvals(v) = table{2,5};
end
truestat = zeros(bsize);
truestat(implicitmask) = fvals;
tfcestat = transform(truestat,H,E,C,dh);
tfcestat = tfcestat(implicitmask);

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute labels
    roccimgs = occimgs;
    for s = 1:nsub
        relabeling = randsample(levels,levels);
        roccimgs(s,relabeling,:) = roccimgs(s,:,:);
    end
    
    % calculate permutation statistic
    rfvals = NaN(nvox,1);
    for v = 1:nvox
        [~,table] = anova_rm(roccimgs(:,:,v),'off');
        rfvals(v) = table{2,5};
    end
    rbrain = zeros(bsize);
    rbrain(implicitmask) = rfvals;
    rbrain = transform(rbrain,H,E,C,dh);
    rstats = rbrain(implicitmask);
    
    % compare maxima to true statistic and increment as appropriate
    curexceeds = max(rstats) >= tfcestat;
    exceedances = exceedances + curexceeds;
end

% create corrected p-value image
corrected = exceedances./nperm;
pcorr = ones(bsize);
pcorr(implicitmask) = corrected;

end

