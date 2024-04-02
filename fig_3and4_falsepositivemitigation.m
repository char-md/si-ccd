close all; clear; clc

% This code is based on the Matlab script files included as part of the 
% GOTCHA CCD challenge problem dataset
%
% The original version of the scripts were written by 
% Steven Scarborough and LeRoy Gorham (AFRL/RYAP)           %
% Email:  steven.scarborough@wpafb.af.mil / leroy.gorham@wpafb.af.mil

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex');


%% Problem Parameters
% Define the path to the base directory of the dataset
datadir = '~/Downloads/ChangeDetectionDataset_GOTCHA/CCD-CP-XPol/HH/';

% Define input pass names here
pass1 = 'FP0120';             % Which pass is the reference image
pass2 = 'FP0124';             % Which pass is the mission image

% Define change detection parameters here
CCDwinSize = 03;              % Window size for CCD (see help file for CCDbasic.m)

% Define image parameters here
minRow = 1300;                % Minimum row value
maxRow = 3300-1;                % Maximum row value
minCol = 1200;                % Minimum column value
maxCol = 3200-1;                % Maximum column value

dyn_range = 60;               % dB of dynamic range to display on SAR image outputs


%% Read in and display data 
% Determine the file names of the input files
im1filename = sprintf('%s%s/c00007a283p50.mat',datadir,pass1);
im2filename = sprintf('%s%s/c00007a283p50.mat',datadir,pass2);

% Load in the reference image
data1 = load(im1filename);

% Crop the image to requested size
im1 = data1.SARdataOut(minRow:maxRow,minCol:maxCol);

% Display the reference image
figure(1)
imagesc([],[],20*log10(abs(im1)./max(abs(im1(:)))),[-dyn_range 0]);
axis image
title('Reference Image');
colorbar

% Load in the mission image
data2 = load(im2filename);

% Crop the image to requested size
im2 = data2.SARdataOut(minRow:maxRow,minCol:maxCol);

% Display the mission image
figure(2)
imagesc([],[],20*log10(abs(im2)./max(abs(im2(:)))),[-dyn_range 0]);
axis image
title('Mission Image');
colorbar


%% Compute change map
% Perform CCD using the MLE method
CCDimage = CCDbasic(im1,im2,CCDwinSize,'MLE',false);
% the last argument determines whether to implement low-RCS masking

% Display the CCD results
figure;
imagesc(CCDimage,[0 1]);
axis image
title('Coherent Change Detection (MLE)');
colormap gray
colorbar


%% Properties of the different regions via the change maps

% Break up change map into regions (square patches)
reg_size = 100;         % size of patch


% Pick three regions - one for the foliage, one for the running track, and
% one for a vehicle insertion/deletion
foliage = CCDimage(101:100+reg_size,101:100+reg_size);
run_track = CCDimage(801:800+reg_size,1751:1751+reg_size);
veh = CCDimage(1001:1000+reg_size,1051:1051+reg_size);


% Plot these regions
figure; 
subplot(2,3,1); imshow(foliage)
axis image
colorbar
title 'Change Map - Foliage'

subplot(2,3,2); imshow(run_track)
axis image
colorbar
title 'Change Map - Running Track'

subplot(2,3,3); imshow(veh)
axis image
colorbar
title 'Change Map - Vehicle Change'




% Analysis (connected components) of the foliage region
% Analysis is done via image connectivity - this uses pixel adjacency to
% determine how many distinct "regions" are there in the image
conn_f = bwconncomp(imbinarize(foliage, 'adaptive', 'ForegroundPolarity', 'bright')); 
                                                     % Compute connectivity
conn_f_lbl = label2rgb(labelmatrix(conn_f), ...      % Label regions for plotting
                        @parula, 'w', 'noshuffle'); 

% Plot these connected regions
subplot(2,3,4); imshow(conn_f_lbl)
axis image
title 'Connected Components - Foliage'


% Analysis of the running track region
conn_r = bwconncomp(imbinarize(run_track, 'adaptive', 'ForegroundPolarity', 'bright')); 
                                                     % Compute connectivity
conn_r_lbl = label2rgb(labelmatrix(conn_r), ...      % Label regions for plotting
                        @parula, 'w', 'noshuffle'); 

% Plot these connected regions
subplot(2,3,5); imshow(conn_r_lbl)
axis image
title 'Connected Components - Running Track'



% Analysis of the vehicle region
conn_v = bwconncomp(imbinarize(veh, 'adaptive', 'ForegroundPolarity', 'bright')); 
                                                     % Compute connectivity
conn_v_lbl = label2rgb(labelmatrix(conn_v), ...      % Label regions for plotting
                        @parula, 'w', 'noshuffle'); 


% Plot these connected regions
subplot(2,3,6); imshow(conn_v_lbl)
axis image
title 'Connected Components - Vehicles'
drawnow;


%% Figure out and plot number of connected components

% Initialize connectivity matrix
ncon = zeros(size(im1));

% Offset for speedup
offset = 5;

% Break up change map into regions (square patches)
reg_size = 30;         % size of patch

for ix = reg_size/2+1:offset:size(im1,1)-reg_size/2
    for iy = reg_size/2+1:offset:size(im1,2)-reg_size/2
        % Extract image
        im = CCDimage(ix-reg_size/2:ix+reg_size/2, iy-reg_size/2:iy+reg_size/2);
        
        % Connectivity
        tmp = bwconncomp(imbinarize(im, 'adaptive', 'ForegroundPolarity', 'bright'));
        ncon(ix-reg_size/2:ix+reg_size/2, iy-reg_size/2:iy+reg_size/2) = tmp.NumObjects;
    end
end

% Plot the number of connected components 
figure;
imagesc(ncon)
axis image
colorbar
title 'No. of connected components in a 31x31 neighborhood'


%% Post-process MLE change map using connectivity information

% Connectivity threshold
th = 0.35*max(ncon(:));

% Generate mask using connectivity information
mask = (ncon<=th);

% Generate post-processed change map
CCDimage_pp = CCDimage;
CCDimage_pp(mask==0) = 1;

% Plot revised change map

figure;
imagesc(CCDimage_pp,[0 1]);
axis image
title('Coherent Change Detection (MLE+connectivity post-processing)');
colormap gray
colorbar