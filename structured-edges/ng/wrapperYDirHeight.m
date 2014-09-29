function wrapperYDirHeight(inName, outName)
  % Load the point cloud from somewhere
  z = imread(inName);
  C = cropCamera(getCameraParam('color'));

  % Convert to cm to pass into the computeCues function
  z = double(z)./10;
  D = z;
  
  [y1 y2 y3 angl1 angl2] = yCues(D, C, 1);
  save(outName, 'y1', 'y2', 'y3', 'angl1', 'angl2');
end
