function wrapperNormalCues(inName, outName)
  % Load the image
  z = imread(inName);

  % Convert to cm to pass into the computeCues function
  z = double(z)./10;

  % The camera matrix from somewhere
  C = cropCamera(getCameraParam('color'));
  
  [dt.ng1, dt.ng2, dt.dg] = normalCues(z, C, 1);
  save(outName, '-STRUCT', 'dt');
end
