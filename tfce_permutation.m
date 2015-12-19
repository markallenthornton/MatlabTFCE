function [pcorr] = tfce_permutation(imgs,varargin)
%TFCE_PERMUTATION performs stepdown maximal statistic permutation testing
%   [pcorr] = tfce_permutation(imgs,nperm) corrects a one-sample,
%   one-tailed t-test for multiple comparisons via a stochastic permutation
%   procedure (random sign flipping). Maximal t-statistics from tests on
%   permuted data are compared with t-statistics in the unshuffled original
%   data using a stepdown procedure described by Holmes, Blair, Watson, and
%   Ford (1996) with an implementation due to Westfall and Young (1993).
%   This procedure is analogous to the more commonly known sequentially
%   rejective test due to Holm (1979).
%
%   Arguments:
%   imgs -- a 4D (x,y,z,subject) matrix of images
%   nperm -- number of permutations to perform. More permutations yield
%   more precise correct p-values. Default set to 1000, but 10000 suggested
%   for publication purposes.
%
%   Output:
%   pcorr -- wholebrain map of corrected p-values

% set defaults
nperm = 1000;
if nargin > 1
    nperm = varargin{1};
end

% calculate matrix size
bsize = size(imgs);
nsub = bsize(4);
bsize = bsize(1:3);

% calculate true t-statistic image
truestat = mean(imgs,4)./(std(imgs,0,4)./sqrt(nsub));
implicitmask = ~isnan(truestat);
tfcetrue = tfce_transform(truestat);


% sort p-values for comparison
tvals = tfcetrue(implicitmask);
[stvals,tind] = sort(tvals,1,'descend');
nvox = length(tvals);

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute signs
    relabeling = randsample([-1 1],nsub,'true');
    rimgs = imgs;
    for s = 1:nsub
        if relabeling(s) == -1;
            rimgs(:,:,:,s) = -imgs(:,:,:,s);
        end
    end
    
    % calculate permutation t-values
    rstat = mean(rimgs,4)./(std(rimgs,0,4)./sqrt(nsub));
    rtfce = tfce_transform(rstat);
    rtvals = rtfce(implicitmask);
    rtvals = rtvals(tind);
    
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

