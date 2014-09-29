
addpath /home/fowlkes/gpb_owt/grouping/lib/
outFile = '/home/fowlkes/gpb_owt/grouping/data/101087_gPb.mat_pbs.mat';
simFile = '/home/fowlkes/gpb_owt/grouping/data/101087_gPb.mat_pbs.mat_sim.tmp';
imgFile = '/home/fowlkes/gpb_owt/grouping/data/101087_gPb.mat_pbs.mat_img.jpg';
latFile = '/home/fowlkes/gpb_owt/grouping/data/101087_gPb.mat_pbs.mat_lat.tmp';
load /home/fowlkes/gpb_owt/grouping/data/mPb.mat  %load in mPb data.

%
% old code
%
l{1} = zeros(size(mPb, 1) + 1, size(mPb, 2));
l{1}(2:end, :) = mPb;
l{2} = zeros(size(mPb, 1), size(mPb, 2) + 1);
l{2}(:, 2:end) = mPb;
latFile = strcat(outFile, '_lat.tmp');
write_array(latFile, l);

system(sprintf('/home/fowlkes/gpb_owt/grouping/lib/segment -image %s -smatrixfile %s -readlattice true -latticefile %s -dthresh 5', imgFile, simFile, latFile));
Wold = read_smatrix(simFile);


%
% new code
%

% mex buildW.cpp -Iutil smatrix.cc ic.cc affinity.cc util/libutil.a

[x,y,z] = buildW(l{1},l{2});
Wnew = sparse(x,y,z);

sum(sum( (Wold-Wnew).^2 ))


  
