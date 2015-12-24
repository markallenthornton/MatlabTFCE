function [tfced] = matlab_tfce_transform_twotailed(img,H,E,C,dh)
%MATLAB_TFCE_TRANSFORM_TWOTAILED convenience wrapper for tfce_transform
%   Performs threshold free cluster enhancement of both positive and
%   negative sides of an image.
%   -- img the 3D image to be transformed
%   -- H height exponent
%   -- E extent exponent
%   -- C connectivity
%   -- dh step size for cluster formation

% call tfce on both sides of the image
tfced_pos = matlab_tfce_transform(img,H,E,C,dh);
tfced_neg = matlab_tfce_transform(-img,H,E,C,dh);
tfced = tfced_pos-tfced_neg;

end

