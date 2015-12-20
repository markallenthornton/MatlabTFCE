function [pcorr_pos,pcorr_neg] = tfce_permutation_independent_twotailed(imgs1,imgs2,varargin)
%TFCE_PERMUTATION_INDEPENDENT_TWOTAILED two-tailed verison of independent
%samples means comparison (see tfce_permutation_independent.m).

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
[stvals,tind] = sort(abs(tvals),1,'descend');
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

