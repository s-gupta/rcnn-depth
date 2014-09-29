### Installation Instructions: ###
0. Checkout eccv14-code, utils, rgbdutils, nyu-hooks
  ```
  #!sh
  git clone git@bitbucket.org:saurabhgupta/eccv14-code.git
  git clone git@bitbucket.org:saurabhgupta/rgbdutils.git eccv14-code/rgbdutils
  git clone git@bitbucket.org:saurabhgupta/utils.git eccv14-code/utils
  ```
0. Checkout caffe-code 
  ```
  #!sh
  git clone https://github.com/BVLC/caffe.git eccv14-code/caffe
  cd eccv14-code/caffe
  git checkout e5cc609138a0bc4ce5177a67cf84952756d11b38
  cd ../../
  ```
0. Get the data (color image, depth images, rawdepth images, splits, ground truth, tasks)
  ```
  #!sh
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-data.tgz
  tar -xf eccv14-data.tgz
  ```
0. Get the caffe models
  ```
  #!sh
  cd eccv14-code
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-caffe-data.tgz
  tar -xf eccv14-caffe-data.tgz 
  cd ..
  ```
0. Set the following paths:
  0. in eccv14-code/get_paths.m : p.root_dir (directory where all intermediate results will be stored).

### Building ###
* Build caffe
* Build MCG, RCNN
```
#!matlab
mcg_build();
rcnn_build();
```

### Edges, UCMs and region proposals ###
* Run the following on the MATLAB prompt. 
```
#!matlab
jobName = 'edges_to_ucms'; script_regions;
jobName = 'benchmark_multi_ucm'; script_regions;
jobName = 'pareto'; script_regions;
jobName = 'cache-mcg'; script_regions;
jobName = 'rank_training'; script_regions;
jobName = 'region-detect'; script_regions;
```

### Detector finetuning and training ###
* Run the following on the MATLAB prompt. At the end of this, you will get 2 finetuning commands that you need to use with caffe for finetuning.
```
#!matlab
jobName = 'save_color'; script_detection;
jobName = 'save_hha'; script_detection;
jobName = 'write_window_file'; script_detection;
```
* Use the finetuning commands to finetune the CNN. Use the following to extract features and train the RCNN model.
```
#!matlab
jobName = 'hha_cache_features'; script_detection;
jobName = 'color_cache_features'; script_detection;
res = rcnn_all('task-guptaetal', 'train', 'val');
```
