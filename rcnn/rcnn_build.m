function rcnn_build()
cd('rcnn');
% Compile liblinear
if ~exist('liblinear_train')
  fprintf('Compiling liblinear version 1.94\n');
  fprintf('Source code page:\n');
  fprintf('   http://www.csie.ntu.edu.tw/~cjlin/liblinear/\n');
  mex -outdir . ...
      CFLAGS="\$CFLAGS -xc++ -O3 -fPIC" -largeArrayDims ...
      external/liblinear-1.94/matlab/train.c ...
      external/liblinear-1.94/matlab/linear_model_matlab.c ...
      external/liblinear-1.94/linear.cpp ...
      external/liblinear-1.94/tron.cpp ...
      "external/liblinear-1.94/blas/*.c" ...
      -output liblinear_train;
  cd('..');
end
