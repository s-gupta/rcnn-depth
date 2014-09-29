#### Installation Instructions ####
0. Create directory, checkout eccv14-code, utils, rgbdutils, nyu-hooks

  ```shell
  mkdir rcnn-depth && cd rcnn-depth
  git clone git@github.com:s-gupta/rcnn-depth.git eccv14-code
  git clone git@github.com:s-gupta/rgbdutils.git eccv14-code/rgbdutils
  git clone git@github.com:s-gupta/utils.git eccv14-code/utils
  git clone git@github.com:s-gupta/nyu-hooks.git eccv14-code/nyu-hooks
  ```

0. Checkout caffe-code 

  ```shell
  git clone https://github.com/BVLC/caffe.git eccv14-code/caffe
  cd eccv14-code/caffe
  git checkout e5cc609138a0bc4ce5177a67cf84952756d11b38
  cd ../../
    ```
  
0. Get the data (color image, depth images, rawdepth images, splits, ground truth, tasks), and external model data (Caffe trained Imagenet model, structured forests BSDS model).

  ```shell
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-data.tgz
  tar -xf eccv14-data.tgz
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-external-data.tgz
  tar -xf eccv14-external-data.tgz
  ```

0. Get precomputed models.
 
  ```
  cd eccv14-code
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-models.tgz
  tar -xf eccv14-models.tgz 
  cd ..
  ```

### Building ###
0. Build caffe (Adjust paths for CUDA / MATLAB in Makefile.config.example and copy to Makefile.config)

  ```shell
  cd eccv14-code/caffe
  make -j 16
  make -j 16 matcaffe
  cd ../..
  ```
  
0. Build imagestack (Adjust paths in eccv14-code/rgbdutils/imagestack/Makefile).
  ```shell
  cd eccv14-code/rgbdutils/imagestack/
  make all -j 16
  cd ..
  ```
  
0. Build toolboxes, MCG, RCNN.

  ```matlab
  mcg_build();
  rcnn_build();
  structured_edges_build();
  
  ```

### Inference ###
#### Contour Detection, UCMs, Region Proposals and Detection on a new image ####
  ```matlab
  %%
  demo();
  
  %% See run_all.m
  [E, ucm2, candidates, detection_scores_no_nms, cls] = ...
    run_all(color_image, depth_image, rawdepth_image, camera_matrix, []);
  ```
  
### Training ###
#### Contour Detection ####
0. Run the following in MATLAB

  ```matlab
  jobName = 'compute_edge_cues'; script_edges;
  jobName = 'train_edge_model'; script_edges;
  jobName = 'test_edge_model'; script_edges;
  ```

#### UCMs and Region Proposals ####
0. Run the following in MATLAB

  ```matlab
  jobName = 'edges_to_ucms'; script_regions;
  jobName = 'benchmark_multi_ucm'; script_regions;
  jobName = 'pareto'; script_regions;
  jobName = 'cache-mcg'; script_regions;
  jobName = 'rank_training'; script_regions;
  jobName = 'region-detect'; script_regions;
  ```

#### Object detectors: finetuning, training ####
0. Run the following in MATLAB. At the end of this, you will get 2 finetuning commands that you need to use with caffe for finetuning.

  ```matlab
  jobName = 'save_color'; script_detection;
  jobName = 'save_hha'; script_detection;
  jobName = 'write_window_file'; script_detection;
  ```

0. Use the finetuning commands that print out to finetune the CNN. Use the following to extract features and train the RCNN model.

  ```matlab
  jobName = 'hha_cache_features'; script_detection;
  jobName = 'color_cache_features'; script_detection;
  res = rcnn_all('task-guptaetal', 'train', 'val');
  ```
