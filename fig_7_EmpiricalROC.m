% Script to generate simulated change detection data

close all; clear; clc

% for pretty pictures
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultLegendInterpreter','latex');


% Set seed for repeatability
rng(31415);


%% Parameters
% Scene/image size
nx = 200;
ny = 205;


% Grid indices
[xj, yj] = ndgrid((0:nx-1).', (0:ny-1).');


% We will have two fixed/big objects and three small/moved objects; these
% small objects will be our surrogate for scene changes

% Position of fixed (big) objects
% These are rectangular; the co-ordinates below correspond to the centers
% of the objects
nbx = [40; 140];        % x-coordinate of the big objects
nby = [75; 150];        % y-coordinate of the big objects


% Position of moved (smallg) objects
% These are centers of the objects
nsx = [80; 160; 130; 100; 50; 160];     % x-coordinate of the small objects
nsy = [50; 100; 60; 70; 160; 40];       % y-coordinate of the small objects


% Locations where ROC is computed
% The first two are locations where there is a change;
% The last two are locations of no change
roc_x = [50; 40; 75; 120];
roc_y = [80; 160; 40; 80];


% No. of ROC trials
ntrials = 1000;

% No. of thresholds in ROC curve
nthr = 1000;

% Noise variance
sig = sqrt(1e-1);


%% Generate time-1 image
im1 = ones(nx,ny);

im1true = im1;          % For calculating image SNR

% Background noise
% This is circularly symmetric complex normnal
im1 = im1 + raylrnd(sig,nx,ny).*exp(2i*pi*rand(nx,ny));



% Add the big objects
im1(nbx(1)-10:nbx(1)+10, nby(1)-50:nby(1)+50) = 2 + ...
                        raylrnd(1e-2,21,101).*exp(2i*pi*rand(21,101));
im1(nbx(2)-40:nbx(2)+40, nby(2)-10:nby(2)+10) = 3 + ...
                        raylrnd(1e-2,81,21).*exp(2i*pi*rand(81,21));

% For SNR computation; add objects to true image
im1true(nbx(1)-10:nbx(1)+10, nby(1)-50:nby(1)+50) = 2;
im1true(nbx(2)-40:nbx(2)+40, nby(2)-10:nby(2)+10) = 3;


% Add the small objects
% The first two objects are square
im1(nsx(1)-5:nsx(1)+5, nsy(1)-5:nsy(1)+5) = 2.5 + ...
                        raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
im1(nsx(2)-5:nsx(2)+5, nsy(2)-5:nsy(2)+5) = 3.5 + ...
                        raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
% The third object is circular
idx = ((xj-nsx(3)).^2 + (yj-nsy(3)).^2 <= 4^2);
im1(idx) = im1(idx) + 4.5;

% For SNR computation; add objects to true image
im1true(nsx(1)-5:nsx(1)+5, nsy(1)-5:nsy(1)+5) = 2.5;
im1true(nsx(2)-5:nsx(2)+5, nsy(2)-5:nsy(2)+5) = 3.5;
im1true(idx) = im1true(idx) + 4.5;



% Here is the SNR
[psnr, snr] = psnr(real(im1), im1true);



% Plot magnitude image
figure; 
imagesc(abs(im1), [0,5]); hold on
axis image
colorbar
title ('Time-1 image (magnitude)')

% Indicate points where ROC is computed
plot3(roc_x(1), roc_y(1), 5, 'rx')
plot3(roc_x(2), roc_y(2), 5, 'rx')
plot3(roc_x(3), roc_y(3), 5, 'rx')
plot3(roc_x(4), roc_y(4), 5, 'rx')


%% Generate time-2 image
im2 = ones(nx,ny);


% Background noise
% This is circularly symmetric complex normnal
im2 = im2 + raylrnd(sig,nx,ny).*exp(2i*pi*rand(nx,ny));


% Add the big objects
im2(nbx(1)-10:nbx(1)+10, nby(1)-50:nby(1)+50) = 2 + ...
                (sig/1e1)*raylrnd(1e-2,21,101).*exp(2i*pi*rand(21,101));
im2(nbx(2)-40:nbx(2)+40, nby(2)-10:nby(2)+10) = 3 + ...
                (sig/1e1)*raylrnd(1e-2,81,21).*exp(2i*pi*rand(81,21));


% Add the small objects
% The first two objects are square
im2(nsx(4)-5:nsx(4)+5, nsy(4)-5:nsy(4)+5) = 2.5 + ...
                    (sig/1e1)*raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
im2(nsx(5)-5:nsx(5)+5, nsy(5)-5:nsy(5)+5) = 3.5 + ...
                    (sig/1e1)*raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
% The third object is circular
idx = ((xj-nsx(6)).^2 + (yj-nsy(6)).^2 <= 4^2);
im2(idx) = im2(idx) + 4.5;



% Plot magnitude image
figure; 
imagesc(abs(im2), [0,5]); hold on
axis image
colorbar
title ('Time-2 image (magnitude)')

% Indicate points where ROC is computed
plot3(roc_x(1), roc_y(1), 5, 'rx')
plot3(roc_x(2), roc_y(2), 5, 'rx')
plot3(roc_x(3), roc_y(3), 5, 'rx')
plot3(roc_x(4), roc_y(4), 5, 'rx')


%% Simple change map

% Perform CCD using the MLE method
CCDwinSize = 3;     % Neighborhood size
CCDimage = CCDbasic(im1,im2,CCDwinSize,'MLE',true);
% the last argument determines whether to implement low-RCS masking


% Display the CCD results
figure;
imagesc(CCDimage,[0 1]); hold on
axis image
title('Coherent Change Detection (MLE)');
colormap gray
colorbar


% Indicate points where ROC is computed
plot3(roc_x(1), roc_y(1), 5, 'rx')
plot3(roc_x(2), roc_y(2), 5, 'rx')
plot3(roc_x(3), roc_y(3), 5, 'rx')
plot3(roc_x(4), roc_y(4), 5, 'rx')


%% Change map with edge information

% Filtering and edge detection
im1_filt = specklefilt(abs(im1),DegreeOfSmoothing=0.6,NumIterations=50);
im1_edge = sqrt(ComputeEdgesImage(im1_filt,'xdim',CCDwinSize).^2 + ComputeEdgesImage(im1_filt,'ydim',CCDwinSize).^2 );


% Display edge map
figure; 
imagesc(im1_edge); 
axis image
title('Edge map - time-1 image')
colormap gray
colorbar


im2_filt = specklefilt(abs(im2),DegreeOfSmoothing=0.6,NumIterations=50);
im2_edge = sqrt( ComputeEdgesImage(im2_filt,'xdim',CCDwinSize).^2 + ComputeEdgesImage(im2_filt,'ydim',CCDwinSize).^2 );



% Plot a cross-section of the function magnitudes and edge maps
figure
plot(abs(im1(80,:))); hold on
plot(im1_edge(80,:), ':')
plot(abs(im2(80,:)), '--')
plot(im2_edge(80,:), '-.')
xlabel('$x$-coordinate index $j$'); ylabel('$|f(x_j, y_0)|$')
xlim([1 nx]); grid on
legend('$f$ (time-1)', '$g$ (time-2)', 'Edges in $f$', 'Edges in $g$')
title('Cross-section of function magnitude and edge map')




% Compute change map
CCDimage = CCD_mpol(cat(3,im1,im1_edge), cat(3,im2,im2_edge), CCDwinSize);


% Display the CCD results
figure;
imagesc(CCDimage,[0 1]); hold on
axis image
title('Coherent Change Detection (w. edge info.)');
colormap gray
colorbar


% Indicate points where ROC is computed
plot3(roc_x(1), roc_y(1), 5, 'rx')
plot3(roc_x(2), roc_y(2), 5, 'rx')
plot3(roc_x(3), roc_y(3), 5, 'rx')
plot3(roc_x(4), roc_y(4), 5, 'rx')


%% Compute ROC

% Store results here
pfa = zeros(nthr, 1);
pd = zeros(nthr, 1);


% Generate thresholds
gma = linspace(0,1,nthr).';


for ix = 1:ntrials
    
    % Generate time-1 image
    im1 = generate_time1_image(nx, ny, nbx, nby, nsx, nsy, xj, yj, sig);

    % Generate time-2 image
    im2 = generate_time2_image(nx, ny, nbx, nby, nsx, nsy, xj, yj, sig);

    % Compute change map
    CCDimage = CCDbasic(im1,im2,CCDwinSize,'MLE',true);


    % Check for false alarms
    pfa = pfa + (CCDimage(roc_y(3), roc_x(3))<gma)/(2*ntrials);
    pfa = pfa + (CCDimage(roc_y(4), roc_x(4))<gma)/(2*ntrials);

    % Check for correct detects
    pd = pd + (CCDimage(roc_y(1), roc_x(1))<gma)/(2*ntrials);
    pd = pd + (CCDimage(roc_y(2), roc_x(2))<gma)/(2*ntrials);
    
end


% Plot ROC curve
figure; 
plot(pfa, pd, '-'); hold on
xlabel '$P_{fa}$'; ylabel '$P_d$'
title 'Empirical ROC curve - change detection'
grid on
drawnow



% Now, the ROC with edge information


% Store results here
pfa_e = zeros(nthr, 1);
pd_e = zeros(nthr, 1);



for ix = 1:ntrials
    
    % Generate time-1 image
    im1 = generate_time1_image(nx, ny, nbx, nby, nsx, nsy, xj, yj, sig);

    % Filtering and edge detection
    im1_filt = specklefilt(abs(im1),DegreeOfSmoothing=0.6,NumIterations=50);
    im1_edge = sqrt(ComputeEdgesImage(im1_filt,'xdim',CCDwinSize).^2 + ComputeEdgesImage(im1_filt,'ydim',CCDwinSize).^2 );


    % Generate time-2 image
    im2 = generate_time2_image(nx, ny, nbx, nby, nsx, nsy, xj, yj, sig);

    % Filtering and edge detection
    im2_filt = specklefilt(abs(im2),DegreeOfSmoothing=0.6,NumIterations=50);
    im2_edge = sqrt( ComputeEdgesImage(im2_filt,'xdim',CCDwinSize).^2 + ComputeEdgesImage(im2_filt,'ydim',CCDwinSize).^2 );


    % Compute change map
    CCDimage = CCD_mpol(cat(3,im1,im1_edge), cat(3,im2,im2_edge), CCDwinSize);

    % Check for false alarms
    pfa_e = pfa_e + (CCDimage(roc_y(3), roc_x(3))<gma)/(2*ntrials);
    pfa_e = pfa_e + (CCDimage(roc_y(4), roc_x(4))<gma)/(2*ntrials);

    % Check for correct detects
    pd_e = pd_e + (CCDimage(roc_y(1), roc_x(1))<gma)/(2*ntrials);
    pd_e = pd_e + (CCDimage(roc_y(2), roc_x(2))<gma)/(2*ntrials);
    
end


% Add this ROC to the plot
plot(pfa_e, pd_e, '-'); hold on
legend( 'MLE', 'w. edge info.' );



%% Helper functions

% ... for generating time-1 image
function im1 = generate_time1_image(nx, ny, nbx, nby, nsx, nsy, xj, yj, sig)

im1 = ones(nx,ny);


% Background noise
% This is circularly symmetric complex normnal
im1 = im1 + raylrnd(sig,nx,ny).*exp(2i*pi*rand(nx,ny));


% Add the big objects
im1(nbx(1)-10:nbx(1)+10, nby(1)-50:nby(1)+50) = 2 + ...
                (sig/1e1)*raylrnd(1e-2,21,101).*exp(2i*pi*rand(21,101));
im1(nbx(2)-40:nbx(2)+40, nby(2)-10:nby(2)+10) = 3 + ...
                (sig/1e1)*raylrnd(1e-2,81,21).*exp(2i*pi*rand(81,21));


% Add the small objects
% The first two objects are square
im1(nsx(1)-5:nsx(1)+5, nsy(1)-5:nsy(1)+5) = 2.5 + ...
                    (sig/1e1)*raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
im1(nsx(2)-5:nsx(2)+5, nsy(2)-5:nsy(2)+5) = 3.5 + ...
                    (sig/1e1)*raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
% The third object is circular
idx = ((xj-nsx(3)).^2 + (yj-nsy(3)).^2 <= 4^2);
im1(idx) = im1(idx) + 4.5;


end






% ... for generating time-2 image
function im2 = generate_time2_image(nx, ny, nbx, nby, nsx, nsy, xj, yj, sig)

im2 = ones(nx,ny);


% Background noise
% This is circularly symmetric complex normnal
im2 = im2 + raylrnd(sig,nx,ny).*exp(2i*pi*rand(nx,ny));


% Add the big objects
im2(nbx(1)-10:nbx(1)+10, nby(1)-50:nby(1)+50) = 2 + ...
                (sig/1e1)*raylrnd(1e-2,21,101).*exp(2i*pi*rand(21,101));
im2(nbx(2)-40:nbx(2)+40, nby(2)-10:nby(2)+10) = 3 + ...
                (sig/1e1)*raylrnd(1e-2,81,21).*exp(2i*pi*rand(81,21));


% Add the small objects
% The first two objects are square
im2(nsx(4)-5:nsx(4)+5, nsy(4)-5:nsy(4)+5) = 2.5 + ...
                (sig/1e1)*raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
im2(nsx(5)-5:nsx(5)+5, nsy(5)-5:nsy(5)+5) = 3.5 + ...
                (sig/1e1)*raylrnd(1e-2,11,11).*exp(2i*pi*rand(11,11));
% The third object is circular
idx = ((xj-nsx(6)).^2 + (yj-nsy(6)).^2 <= 4^2);
im2(idx) = im2(idx) + 4.5;



end