function imOut = CCDbasic(imIn1,imIn2,winSize,method,lowRCSmask,lowRCSmaskmethod)
    arguments
        % This defines the expected input parameter data types
        imIn1 (:,:) double 
        imIn2 (:,:) double 
        winSize (1,1) uint8

        % This sets default values for optional parameters
        method (1,1) string = 'MLE'
        lowRCSmask (1,1) logical = false
        lowRCSmaskmethod (1,1) string = 'entropy'
        % lowRCSmaskmethod (1,1) string = 'novak'
    end

% This function is a modified version of the function file made available
% with the GOTCHA CCD dataset
% The original help comments are included below
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs a basic non-coherent change detection (NCCD)  %
% operation.  The inputs are:                                          %
%                                                                      %
% imIn1:  First input image (reference image)                          %
% imIn2:  Second input image (mission image)                           %
% winSize:  Size of window (can be rectangular or square)              %
%     For square window, winSize is a scalar, and the resulting window %
%         has dimensions (2*winSize+1) x (2*winSize+1)                 %
%     For rectangular window, winSize is a 2-vector, and the resulting %
%         window has dimensions (2*winSize(1)+1) x (2*winSize(2)+1)    %
%                                                                      %
% The output is:                                                       %
% imOut:  The CCD coherence map.  The output values are from 0 to 1,   %
%     with 0 indicating no coherence in the two images and 1           %
%     indicating full coherence.                                       %
%                                                                      %
% References:                                                          %
%   Scarborough, S. "A Challenge Problem for SAR Change Detection and  %
%       Data Compression," SPIE Algorithms for Synthetic Aperture      %
%       Radar Imagery XVII, Orlando, FL, April, 2010.                  %
%                                                                      %
%   Novak, L. "Coherent Change Detection for Multi-Polarization SAR,"  %
%       Asilomar Conference on Circuits, Systems, and Computers,       %
%       Pacific Grove, CA, October, 2005.                              %
%                                                                      %
% Contact Information:                                                 %
% Steven Scarborough and LeRoy Gorham (AFRL/RYAP)                      %
% Email:  steven.scarborough@wpafb.af.mil / leroy.gorham@wpafb.af.mil  %                                  %
% Date Released:  8 Apr 2011                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the window function
switch length(winSize)
    case 1
        % winFunc = ones(2*winSize + 1, 2*winSize + 1);
        % winFunc = fspecial('disk', double(winSize));
        winFunc = fspecial('average', double(2*winSize+1));
    case 2
        winFunc = ones(2*winSize(1) + 1, 2*winSize(2) + 1);
    otherwise
        error('winSize must have a length of 1 or 2');
end


% Perform low-RCS (shadow region) masking?
% This method is based on the paper 
% Change Detection Experiments Using Gotcha Public Release SAR Data
% Algorithms for Synthetic Aperture Radar Imagery XX
% Stojanovic, Ivana and Novak, Les,
% vol. 8746, pp. 144 - 153, SPIE, June 2013
% doi: 10.1117/12.2020650

% TODO: include the Entropy masking method

if(lowRCSmask)
    switch lower(lowRCSmaskmethod)
        case 'entropy'
            % Apply entropy filtering
            mask_fnc_im1 = entropyfilt(abs(imIn1),true(2*winSize+1));
            mask_fnc_im2 = entropyfilt(abs(imIn2),true(2*winSize+1));
            % the second argument is the neighborhood size to use
            
            % Reconcile results from the two images
            mask_fnc = min(mask_fnc_im1, mask_fnc_im2);

            % Rescale to [0,1]
            mask_fnc = 1-rescale(mask_fnc);

        case 'novak'
            % This is the equation from Sec. 3 of the Novak manuscript
            mask_fnc = ( conv2(abs(imIn1+imIn2).^2,winFunc,'same') + ...
                        conv2(abs(imIn1-imIn2).^2,winFunc,'same') )/(2*(2*double(winSize)+1));
    end
    
    % Masking threshold
    % TODO: make this an optional parameter input with a default value
    th = 0.15;    

    % Find regions which satisfy the masking condition
    mask_flag = (mask_fnc<=th);

    % Enforce such regions to have unit signal value (so that they don't
    % register in the change map)
    imIn1(mask_flag) = 1;
    imIn2(mask_flag) = 1;
end

% Mean subtraction
imIn1 = imIn1 - conv2(imIn1,winFunc,'same')/numel(winFunc);
imIn2 = imIn2 - conv2(imIn2,winFunc,'same')/numel(winFunc);

% Calculate the numerator of the coherence function
num = abs(conv2(imIn1.*conj(imIn2),winFunc,'same'));

% Depending on the type of CCD, the denominator changes
switch method
    case 'MLE'
        % MLE CCD statistic
        denom = 0.5*(conv2(abs(imIn1).^2,winFunc,'same') + ...
                     conv2(abs(imIn2).^2,winFunc,'same'));
    case 'CCD'
        % "CCD" statistic - with geometric mean normalization
        % Calculate the denominator of the coherence function
        denom = sqrt(conv2(abs(imIn1).^2,winFunc,'same').*...
            conv2(abs(imIn2).^2,winFunc,'same'));
end

% Calculate the coherence function
imOut = num./denom;

return