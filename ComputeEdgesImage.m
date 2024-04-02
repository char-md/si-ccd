function ed = ComputeEdgesImage(im, dim, nbr_size)
% This function computes edges (along the x- or y- dimension) in 2D images

% TODO: Can be made much more efficient with conv2(.) usage!

% Preliminaries
nx = size(im,1);
ny = size(im,2);


% Which dimension are we computing edges along?
switch (dim)
    case 'xdim'
        % Size
        n = ny; 
    case 'ydim'
        % Size
        n = nx;
end


% Here is our edge detection kernel
% We are using divided differences here
% Second order divided differences
D2 = spdiags([-ones(n,1) ones(n,1)], [-1 1], n, n);

% Filtering/smoothing matrix
filtmat = zeros(n);
% Filtering kernel
kern = gausswin(n, 100);

% Normalizing factor?
% TODO: This can be written out analytically
kern = kern/2;

% TODO: This can be simplified by making a circulant (or Toeplitz) matrix
% with shifted rows of the smoothing kernel
% Standard basis
e = zeros(n,1);

for ix = 1:n
    % Standard basis
    e(ix) = 1;

    % Apply filtering
    e_smooth = conv(e, kern, 'same');

    % Column of matrix
    filtmat(:,ix) = e_smooth;

    % Reset entry
    e(ix) = 0;
end


% Include smoothing in the edge detection
D = D2*filtmat;
% D = D2;               % No smoothing


% Here is where we perform the actual filtering
% Output
ed = zeros(size(im));

switch (dim)
    case 'xdim'

        for ix = 1:nx
            ed(ix,:) = (D*im(ix,:).').';
        end

    case 'ydim'

        for iy = 1:ny
            ed(:,iy) = D*im(:,iy);
        end
end


% Perform matched filtering to clean up the response

% Grid
x = linspace(-pi,pi,n).';

% Waveform/template to look for
% This si the response to a unit sized jump
jmp_wvfm = D*(x>=0);

% Extract the neighborhood in the vicinity of jump
[~, cloc] = min(abs(x-0));              % Index of the center/zero location
jmp_wvfm(1:cloc-nbr_size) = 0;
jmp_wvfm(cloc+nbr_size:end) = 0;
norm_fac = norm(jmp_wvfm)^2;

% Perform matched filtering
% For now, set covariance matrix to identity
% TODO: use noise variance and conc. factor to set the covariance matrix

switch (dim)
    case 'xdim'

        for ix = 1:nx
            ed(ix,:) = (conv(ed(ix,:).',jmp_wvfm,'same')/norm_fac).';
        end

    case 'ydim'

        for iy = 1:ny
            ed(:,iy) = conv(ed(:,iy),jmp_wvfm,'same')/norm_fac;
        end
end



end