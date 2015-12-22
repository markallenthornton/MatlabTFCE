function [tfced] = matlab_tfce_transform_twotailed(img,H,E,C,ndh)
%TFCE_TRANSFORM_TWOTAILED convenience wrapper for tfce_transform
%   Performs threshold free cluster enhancement of both positive and
%   negative sides of an image.
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- ndh step number for cluster formation

% call tfce on both sides of the image
tfced = tfce_transform(img,H,E,C,ndh);
tfced = tfce_transform(-tfced,H,E,C,ndh);
tfced = -tfced;

end

