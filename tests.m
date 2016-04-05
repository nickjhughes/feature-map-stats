
% A collection of tests and examples.

%% Image alignment

% Settings
base_im = imread(fullfile(matlabroot, '\toolbox\images\imdata\cameraman.tif'));
n = 10;

% Create mis-aligned stack of images
true_offsets = zeros(n, 2);
images = nan(size(base_im,1), size(base_im,2), n);
images(:,:,1) = base_im;
for j = 2:n
    true_offsets(j,:) = randi([-10, 10], 1, 2);
    images(:,:,j) = imtranslate(base_im, true_offsets(j,:));
end

% Align images to first image in stack
[aligned, offsets] = align_images(images, [], 'full', 'optim');
assert(all(all(-offsets == true_offsets)));

%% Extended spatial decorrelation

% Settings
px = 64;
py = 64;
N = 3;

% Create and mix three sources
[x, y] = meshgrid(1:px, 1:py);
sources = nan(py, px, N);
sources(:,:,1) = sin(2*pi/(px/3)*x).*sin(2*pi/(py/3)*y);
sources(:,:,2) = (sin(2*pi/(px/1.5)*x).*sin(2*pi/(py/1.5)*y)).^3;
sources(:,:,3) = cos(2*pi/(2*px)*x)+cos(2*pi/(2*py)*y);
A = [0.39, -0.56, 0.78; 0.08, 0.44, 0.57; -0.64, -0.95, -0.82];
mixture = reshape((A*reshape(sources, px*py, N).').', [py, px, 3]);

% Separate sources with PCA
x = reshape(mixture, px*py, N);
x = x - repmat(mean(x,1), size(x,1), 1);
[~, score] = pca(x);
pca_result = reshape(score, [py, px, N]);

% Separate sources with ESD
dr = [5, 5];
esd_result = esd(mixture, dr);

% Compare results
figure;
for j = 1:3
    subplot(3,4,4*(j-1)+1);
    imagesc(sources(:,:,j));
end
for j = 1:3
    subplot(3,4,4*(j-1)+2);
    imagesc(mixture(:,:,j));
end
for j = 1:3
    subplot(3,4,4*(j-1)+3);
    imagesc(pca_result(:,:,j));
end
for j = 1:3
    subplot(3,4,4*(j-1)+4);
    imagesc(esd_result(:,:,j));
end
subplot(3,4,1);
title('Sources');
subplot(3,4,2);
title('Mixtures');
subplot(3,4,3);
title('PCA');
subplot(3,4,4);
title('ESD');

%% Orientation distribution

% Settings
n = 128;
sigma = 6;
bins = -5:10:185;
display_bins = bins(2:end-1);

% Construct OP-like map
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);

% Calculate orientation preference distribution
op_dist = orientation_hist(op, true(size(op)), bins);

% Plot
figure;
hold on;
plot([0, 180], 1/length(display_bins)*[1, 1], 'k-');
plot(display_bins, op_dist);
hold off;
xlabel('Orientation');
ylabel('Proportion of Map');

%% Crossing angle distribution

% Settings
n = 128;
sigma = 6;
nangles = 8;
bins = [0:10:80, 91];
display_bins = 5:10:85;

% Construct spatially unrelated OP-like and OD-like maps
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);
od = imfilter(randn(n), h);

% Calculate crossing angle distribution
crossing = crossing_angle_dist(od, op, true(size(od)), nangles, bins);

% Plot
% Note: the distribution for unrelated maps is a sine curve
figure;
hold on;
plot(display_bins, sind(display_bins)/sum(sind(display_bins)), 'k-');
plot(display_bins, crossing);
hold off;
xlabel('Intersection Angle');
ylabel('Proportion of Crossings');

%% Pinwheel-OD relationship

% Settings
n = 128;
sigma = 6;

% Construct spatially unrelated OP-like and OD-like maps
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);
od = imfilter(randn(n), h);

% Calculate pinwheel-OD relationship
pinwod_dist = pinwod(op, od);

% Plot
figure;
bar(1:5, pinwod_dist);
xlabel('OD Bin');
ylabel('Proportion of Pinwheels');

%% Orientation selectivity

% Settings
n = 128;
sigma = 6;

% Construct spatially unrelated OP-like and OD-like maps
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);
od = imfilter(randn(n), h);

% Calculate orientation selectivity in OD bins
selectivity = binned_selectivity(op, od);

% Plot
figure;
plot(1:5, selectivity);
xlabel('OD Bin');
ylabel('Orientation Selectivity');

%% Pinwheel locations

% Settings
n = 128;
sigma = 6;

% Construct OP-like map
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);

% Find pinwheels
pinw = locate_pinwheels(op);

% Display map and pinwheel locations
figure;
hold on;
imagesc(angle(op)/2);
for j = 1:size(pinw, 1)
    plot(pinw(j,1), pinw(j,2), 'k.', 'MarkerSize', 30);
end
hold off;
colormap hsv;
axis equal;
axis tight;

%% Pinwheel density

% Settings
n = 128;
sigma = 6;

% Construct OP-like map
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);

% Calculate pinwheel density
[~, density] = pinw_density(op);
fprintf('Pinwheel density: %.2f pinwheels per square map wavelength\n', density);

%% Plotting OP contours

% Settings
n = 128;
sigma = 6;

% Construct OP-like map
h = fspecial('Gaussian', ceil(10*sigma), sigma);
op = imfilter(randn(n), h) + 1i*imfilter(randn(n), h);

% Plot contours
angles = [-90, -67.5, -45, -22.5, 0, 22.5, 45, 67.5];
c = op_contours(op, angles);
figure;
hold on;
cmap = hsv(length(angles));
for j = 1:length(c)
    plot(c{j}(1,:), c{j}(2,:), 'Color', cmap(j,:));
end
hold off;
axis equal;
axis tight;
