function [pcorr] = tfce_permutation_independent(imgs1,imgs2,varargin)
%TFCE_PERMUTATION_INDEPENDENT tests means difference, independent samples
%   [pcorr] = tfce_permutation_independent(imgs1,imgs2,nperm) performs a
%   comparison of the voxelwise image means of two independent groups (as
%   in an independent t-test). This version performs 1-tailed tests. The 
%   alternative hypothesis is imgs1>imgs2. Maximal statistics from tests on
%   permuted data are compared with statistics in the original
%   data using a stepdown procedure described by Holmes, Blair, Watson, and
%   Ford (1996) with an implementation due to Westfall and Young (1993).
%   This procedure is analogous to the more commonly known sequentially
%   rejective test due to Holm (1979). Note that independent sample tests
%   require homogeneity of variance.
%
%   Arguments:
%   imgs1 -- images from group 1 with dimensions x,y,z,nsubject1
%   imgs1 -- images from group 2 with dimensions x,y,z,nsubject2
%
%   Output:
%   pcorr -- image of p-values, corrected for multiple comparisons

% set defaults
nperm = 1000;
if nargin > 1
    nperm = varargin{1};
end

% calculate matrix size
bsize = size(imgs1);
nsub1 = bsize(4);
bsize = size(imgs2);
nsub2 = bsize(4);
bsize = bsize(1:3);
nsub = nsub1+nsub2;

% calculate true mean image
truestat = mean(imgs1,4)-mean(imgs2,4);
implicitmask = ~isnan(truestat);

% sort p-values for comparison
tvals = truestat(implicitmask);
[stvals,tind] = sort(tvals,1,'descend');
nvox = length(tvals);

% extract occupied voxels for permutation test
imgs = cat(4,imgs1,imgs2);
glabs = [ones(nsub1,1);ones(nsub2,1)*2];
occimgs = NaN(nvox,nsub);
for s = 1:nsub
    curimg = imgs(:,:,:,s);
    occimgs(:,s) = curimg(implicitmask);
end

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute signs
    relabeling = randsample(glabs,nsub);
    rimgs1 = occimgs(:,relabeling==1);
    rimgs2 = occimgs(:,relabeling==2);
    
    % calculate permutation means
    rstats = mean(rimgs1,2)-mean(rimgs2,2);
    rtvals = rstats(tind);

    % calculate maxima
    maxima = zeros(nvox,1);
    maxima(end) = rtvals(end);
    for v = fliplr(1:(nvox-1))
        maxima(v) = max(rtvals(v),maxima(v+1));
    end
    
    % compare maxima to t-values and increment as appropriate
    curexceeds = maxima >= stvals;
    failed = find(curexceeds);
    curexceeds(failed:end) = 1;
    exceedances = exceedances + curexceeds;
end

% create corrected p-value image
corrected = NaN(nvox,1);
corrected(tind) = exceedances./nperm;
pcorr = ones(bsize);
pcorr(implicitmask) = corrected;

end

