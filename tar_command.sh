cd ../
tar -cf eccv14-external-data.tgz eccv14-code/caffe-data eccv14-code/structured-edges/models/forest/ 

cd ../
tar -cf eccv14-models.tgz \ 
  eccv14-cachedir/release/contours/models/forest \
  eccv14-cachedir/release/regions/classifiers \
  eccv14-cachedir/release/regions/features \
  eccv14-cachedir/release/regions/pareto \
  eccv14-cachedir/release/detection/finetuning/protodir \
  eccv14-cachedir/release/detection/finetuning/windowfile \
  eccv14-cachedir/release/detection/finetuning/snapshot/*30000.caffemodel \
  eccv14-cachedir/release/detection/detector/rgb_hha_30000_trainval/rcnn_model.mat \
  eccv14-cachedir/release/detection/detector/rgb_hha_30000_trainval/pr-curves
