function structured_edges_build()
  cd('structured-edges');
  OMPPARAMS = ['-DUSEOMP CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"'];
  s = sprintf('mex private/edgesDetectMex.cpp -outdir private %s', OMPPARAMS); eval(s);
  s = sprintf('mex private/edgesDetectMexNormal.cpp -outdir private %s', OMPPARAMS); eval(s);
  s = sprintf('mex private/edgesNmsMex.cpp -outdir private %s', OMPPARAMS); eval(s);
  s = sprintf('mex private/spDetectMex.cpp -outdir private %s', OMPPARAMS); eval(s);
  s = sprintf('mex private/edgeBoxesMex.cpp -outdir private'); eval(s);
  cd('..');
end
