function [tfced] = matlab_tfce_transform_twotailed(img,varargin)
%TFCE_TRANSFORM_TWOTAILED convenience wrapper for tfce_transform
%   Performs threshold free cluster enhancement of both positive and
%   negative sides of an image.
%   -- img the 3D image to be transformed
%   -- H height exponent, default = 2
%   -- E extent exponent, default = 0.5
%   -- C connectivity, default = 6 (6 = surface, 18 = edge, 26 = corner)
%   -- ndh step number for cluster formation, default = 100
%   More steps will be more precise but will require more time and memory.
%   The H & E default parameter settings match FSL's randomise/fslmaths.
%   The C default setting matches FSL's ramdomise default setting. To
%   match SPM's default cluster forming, use 18 instead. The transformed
%   image is returend as 'tfced'.

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
tfced = -tfced;

end

