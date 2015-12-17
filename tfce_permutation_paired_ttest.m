function [varargout] = tfce_permutation_paired_ttest(imgs1,imgs2,varargin)
% TFCE_PERMUTATION_PAIRED_TTEST convenience wrapper for performing paired
% t-tests on imaging data via calculating image difference and executing
% one-sample t-tests on said difference.
% imgs1 -- measure 1 images
% imgs2 -- measure 2 images
% nperm -- number of permutations (default 1000)
% tails -- 1 or 2 tailed ttest (default 1)
% H -- height exponent, default = 2
% E -- extent exponent, default = 0.5
% C -- connectivity, default = 6 (6 = surface, 18 = edge, 26 = corner)
% ndh -- step number for cluster formation, default = 100
%
% For 1-tailed tests, returned pcorr, for two-tailed tests, returns
% pcorr_pos and pcorr_neg.

% set defaults
nperm = 1000;
tails = 1;
H = 2;
E = .5;
C = 6;
ndh = 100;

% overriding defaults as necessary
if nargin > 2
    nperm = varargin{1};
end
if nargin > 3
    tails = varargin{2};
end
if nargin > 4
    H = varargin{3};
end
if nargin > 5
    E = varargin{4};
end
if nargin > 6
    C = varargin{5};
end
if nargin > 7
    ndh = varargin{6};
end

% calculate image difference
imgsdiff = imgs1-imgs2;

% perform tfce and permutation (note that image should be transformed in
% two-tailed manner even for one-tailed test)
tfced = tfce_transform_twotailed(imgsdiff,H,E,C,ndh);

% perform appropriate permutation test
if tails == 1
    pcorr = tfce_permutation(tfced,nperm);
    varargout{1} = pcorr;
else
    [pcorr_pos,pcorr_neg] = tfce_permutation_twotailed(tfced,nperm);
    varargout{1} = pcorr_pos;
    varargout{1} = pcorr_neg;
end

end