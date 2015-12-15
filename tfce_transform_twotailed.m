function [tfced] = tfce_transform_twotailed(img,varargin)
%TFCE_TRANSFORM_TWOTAILED convenience wrapper for tfce_transform
%   Performs treshold free cluster enhancement of both positive and
%   negative sides of an image (e.g. in preparation for a two-tailed test).
%   For arguments and details, see tfce_transform.m

% setting defaults
H = 2;
E = .5;
C = 6;
ndh = 100;

% overriding defaults as necessary
if nargin > 1
    H = varargin{1};
end
if nargin > 2
    E = varargin{2};
end
if nargin > 3
    C = varargin{3};
end
if nargin > 4
    ndh = varargin{4};
end

% call tfce on both sides of the image
tfced = tfce_transform(img,H,E,C,ndh);
tfced = tfce_transform(-tfced,H,E,C,ndh);
tfced = - tfced;

end

