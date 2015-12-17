function [pcorr_pos,pcorr_neg] = tfce_permutation_twotailed(imgs,varargin)
%TFCE_PERMUTATION_TWOTAILED two-tailed version of tfce_permutation.m
% See that file for details. Note that two corrected images are returned -
% one for effects in the positive direction and the other for effects in
% the negative direction.

% set defaults
nperm = 1000;
if nargin > 1
    nperm = varargin{1};
end

% calculate matrix size
bsize = size(imgs);
nsub = bsize(4);
bsize = bsize(1:3);

% calculate trut t-statistic image
truestat = mean(imgs,4)./(std(imgs,0,4)/sqrt(nsub));
implicitmask = ~isnan(truestat);

% sort p-values for comparison
tvals = truestat(implicitmask);
[stvals,tind] = sort(abs(tvals),1,'descend');
nvox = length(tvals);

% extract occupied voxels for permutation test
occimgs = NaN(nvox,nsub);
for s = 1:nsub
    curimg = imgs(:,:,:,s);
    occimgs(:,s) = curimg(implicitmask);
end

% preallocate maxinde

% cycle through permutations
exceedances = zeros(nvox,1);
for p = 1:nperm
    
    % permute signs
    relabeling = randsample([-1 1],nsub,'true');
    roccimgs = occimgs;
    for s = 1:nsub
        if relabeling(s) == -1;
            roccimgs(:,s) = -occimgs(:,s);
        end
    end
    
    % calculate permutation t-values
    rstats = mean(roccimgs,2)./(std(roccimgs,0,2)/sqrt(nsub));
    rtvals = abs(rstats(tind));
    
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

% split into positive and negative effects
pos = truestat>0;
pcorr_pos = pcorr;
pcorr_pos(~pos) = 1;
pcorr_neg = pcorr;
pcorr_neg(pos) = 1;

end

