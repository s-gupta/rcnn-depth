I = imread(fullfile('demo-data', 'images.png'));
D = imread(fullfile('demo-data', 'depth.png'));
RD = imread(fullfile('demo-data', 'rawdepth.png'));
C = cropCamera(getCameraParam('color'));
out_file = fullfile('demo-data', 'output.mat');
run_all(I, D, RD, C, out_file);
