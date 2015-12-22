function [pcorr] = tfce_correlation(imgs,covariate,varargin)
% TFCE_CORRELATION computes threshold free cluster enhancement for a an
% individual difference correlation between a covariate and brain activity.
% The maximal statistic technique is augmented with sequential stepdown to
% maximize power.
%
%   Arguments:
%   imgs -- a 4D (x,y,z,subject) matrix of images. Like tfce_permutation,
%   these should have had tfce_transform applied already.
%   covariate -- a vector of length = number of subjects containing values
%   to be correlated with brain activity
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values. Default set to 1000, but 10000 suggested
%   for publication purposes.
%
%   Output:
%   pcorr -- wholebrain map of corrected p-values

% set defaults
nperm = 1000;
if nargin > 2
    nperm = varargin{1};
end

% calculate matrix size
bsize = size(imgs);
nsub = bsize(4);
covariate = covariate(:);
bsize = bsize(1:3);

% calculate implicit mask
sumimg = sum(imgs,4);
implicitmask = ~isnan(sumimg) & sumimg~=0;
nvox = sum(implicitmask(:));

% extract occupied voxels for permutation test
occimgs = NaN(nsub,nvox);
for s = 1:nsub
    curimg = imgs(:,:,:,s);
    occimgs(s,:) = curimg(implicitmask);
end

% calculate true correlation image
truestat = corr(occimgs,covariate);

% sort p-values for comparison
[scvals,cind] = sort(truestat,1,'descend');

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute covariates
    rcov = randsample(covariate,nsub);
    
    % calculate permutation correlations
    rstats = corr(occimgs,rcov);
    rcvals = rstats(cind);
    
    % calculate maxima
    maxima = zeros(nvox,1);
    maxima(end) = rcvals(end);
    for v = fliplr(1:(nvox-1))
        maxima(v) = max(rcvals(v),maxima(v+1));
    end
    
    % compare maxima to t-values and increment as appropriate
    curexceeds = maxima >= scvals;
    failed = find(curexceeds);
    curexceeds(failed:end) = 1;
    exceedances = exceedances + curexceeds;
end

% create corrected p-value image
corrected = NaN(nvox,1);
corrected(cind) = exceedances./nperm;
pcorr = ones(bsize);
pcorr(implicitmask) = corrected;

end