function [mssim, ssim_map] = ssim_index3d(img1, img2, sw)

%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author was with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University, USA. He is
%currently with Department of Electrical and Computer Engineering,
%University of Waterloo, Canada.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 4, Apr. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim ssim_map] = ssim_index(img1, img2, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================


if (nargin < 2 || nargin > 5)
    mssim = -Inf;
    ssim_map = -Inf;
    return;
end

if (size(img1) ~= size(img2))
    mssim = -Inf;
    ssim_map = -Inf;
    return;
end

s = size(img1);
% sw = [2 2 2];

max_val = max(img1(:));
img1 = img1 / max_val * 255;
img2 = img2 / max_val * 255;

if (nargin == 2)
    sw = [2, 2, 2];
    ind = find(img1 ~= 0);
    %    window = fspecial('gaussian', 11, 1.5);	%
    window = gkernel(1.5, sw); %
    K(1) = 0.01; % default settings
    K(2) = 0.03; %
    L = 255; %
end

if (nargin == 3)
    if ((s(1) < sw(1)) || (s(2) < sw(2)) || (s(3) < sw(3)))
        mssim = -Inf;
        ssim_map = -Inf;
        return
    end
    %    window = fspecial('gaussian', 11, 1.5);
    window = gkernel(1.5, sw); %
    ind = find(img1 ~= 0);
    K(1) = 0.01; % default settings
    K(2) = 0.03;
    L = 255;
    if (length(K) == 2)
        if (K(1) < 0 || K(2) < 0)
            mssim = -Inf;
            ssim_map = -Inf;
            return;
        end
    else
        mssim = -Inf;
        ssim_map = -Inf;
        return;
    end
end


C1 = (K(1) * L)^2;
C2 = (K(2) * L)^2;
window = window / sum(window(:));
img1 = double(img1);
img2 = double(img2);

% mu1   = convn( img1,window, 'same');
% mu2   = convn( img2,window, 'same');

mu1 = convolution3D_FFTdomain(img1, window);
mu2 = convolution3D_FFTdomain(img2, window);

mu1_sq = mu1 .* mu1;
mu2_sq = mu2 .* mu2;
mu1_mu2 = mu1 .* mu2;
% sigma1_sq = convn( img1.*img1,window, 'same') - mu1_sq;
% sigma2_sq = convn( img2.*img2,window, 'same') - mu2_sq;
% sigma12 = convn( img1.*img2,window, 'same') - mu1_mu2;


sigma1_sq = convolution3D_FFTdomain(img1.*img1, window) - mu1_sq;
sigma2_sq = convolution3D_FFTdomain(img2.*img2, window) - mu2_sq;
sigma12 = convolution3D_FFTdomain(img1.*img2, window) - mu1_mu2;


if (C1 > 0 && C2 > 0)
    ssim_map = ((2 * mu1_mu2 + C1) .* (2 * sigma12 + C2)) ./ ((mu1_sq + mu2_sq + C1) .* (sigma1_sq + sigma2_sq + C2));
else
    numerator1 = 2 * mu1_mu2 + C1;
    numerator2 = 2 * sigma12 + C2;
    denominator1 = mu1_sq + mu2_sq + C1;
    denominator2 = sigma1_sq + sigma2_sq + C2;
    ssim_map = ones(size(mu1));
    index = (denominator1 .* denominator2 > 0);
    ssim_map(index) = (numerator1(index) .* numerator2(index)) ./ (denominator1(index) .* denominator2(index));
    index = (denominator1 ~= 0) & (denominator2 == 0);
    ssim_map(index) = numerator1(index) ./ denominator1(index);
end

temp = zeros(size(img1));
temp(ind) = 1;
iind = temp == 0;
ssim_map(iind) = 1;
mssim = mean(ssim_map(ind));

return

function [gaussKernel] = gkernel(sigma, sk)

% Pierrick Coupe - pierrick.coupe@gmail.com
% Brain Imaging Center, Montreal Neurological Institute.
% Mc Gill University
%
% Copyright (C) 2008 Pierrick Coupe

for x = 1:(2 * sk(1) + 1)
    for y = 1:(2 * sk(2) + 1)
        for z = 1:(2 * sk(3) + 1)
            radiusSquared = (x - (sk(1) + 1))^2 + (y - (sk(2) + 1))^2 + (z - (sk(3) + 1))^2;
            gaussKernel(x, y, z) = exp(-radiusSquared/(2 * sigma^2));
        end
    end
end


function [outVol] = convolution3D_FFTdomain(inVol, inKer)

%% convolution3D_FFTdomain - Performs a fast 3D convolution between a volume and a kernel using mutliplication in the Fourrier space.
%
% Syntax:  [outVol] = convolutionInFFTdomain(inVol,inKer,inMsg)
% Inputs:
%       inVol - input volume (real / complex)
%       inKer - input kernel (real / complex)
% Outputs:
%       outVol - output convolved volume (real / complex) - precision of the output format is the same as the input
%       volume. The output volume is the central part of the convolution with same size as inVol.
%       size(outVol)=size(inVol) ('same' option of convn).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% Author: Christopher Coello
% Work address: Preclinical PET/CT Unit
% Email address: s.c.coello@medisin.uio.no
% Website: http://www.med.uio.no/imb/english/services/public/pet/
% 2012/03/13; Last revision: 2012/03/16
% Created with Matlab version: 7.13.0.564 (R2011b)

%%
% Check if both inputs are real numbers
realInput = isreal(inVol) && isreal(inKer);
% Check the precision of the input volume
if isa(inVol, 'double')
    indDb = 1;
else
    indDb = 0;
end
% Size of the input volumes
inVolSize = size(inVol);
inVolSide = max(inVolSize);
inKerSize = size(inKer);
inKerSide = max(size(inKer));
% Fourrier tranform of the volume and inKer.
% extr(1:3)={};
extr = cell(1, 3);
for iDim = 1:3
    inVol = fft(inVol, inVolSide+inKerSide-1, iDim);
    inKer = fft(inKer, inVolSide+inKerSide-1, iDim);
    extr{iDim} = ceil((inKerSize(iDim) - 1)/2) + (1:inVolSize(iDim));
end
% Multiplication of the Fourrier tranforms
conv_FFT = inVol .* inKer;
% Inverse Fourrier Transform of the convolution
for iDim = (1:3)
    conv_FFT = ifft(conv_FFT, [], iDim);
end
% Crop the side of the image in relation to the size of the kernel
convinVol = conv_FFT(extr{:});
% limit the results
if realInput
    convinVol = real(convinVol);
end
% Format the output into single precision if needed
if indDb
    outVol = convinVol;
else
    outVol = single(convinVol);
end