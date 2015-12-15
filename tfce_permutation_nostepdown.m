function [pcorr] = tfce_permutation_nostepdown(imgs,varargin)
%TFCE_PERMUTATION_NOSTEPDOWN a version of tfce_permutation without the
% stepdown functionality.

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
    
    % compare maxima to t-values and increment as appropriate
    exceedances = exceedances + (max(rstats) >= tvals);
end

% create corrected p-value image
corrected = exceedances./nperm;
pcorr = ones(bsize);
pcorr(implicitmask) = corrected;

end

