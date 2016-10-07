#### _For more recent work that is faster and more accurate, please refer to our latest paper on [Cross Modal Distillation for Supervision Transfer](http://arxiv.org/pdf/1507.00448.pdf). Code and pre-trained models are available:  [detection](https://github.com/s-gupta/fast-rcnn/tree/distillation), [instance segmentation](https://github.com/s-gupta/fast-rcnn/tree/sds-distillation)._


###Learning Rich Features from RGB-D Images for Object Detection and Segmentation ([pdf](http://www.cs.berkeley.edu/~sgupta/pdf/rcnn-depth.pdf))###
*Saurabh Gupta, Ross Girshick, Pablo Arbeláez and Jitendra Malik*

Presented at European Conference on Computer Vision (ECCV), 2014 

In this paper we study the problem of object detection for RGB-D images using semantically rich image and depth features. We pro- pose a new geocentric embedding for depth images that encodes height above ground and angle with gravity for each pixel in addition to the horizontal disparity. We demonstrate that this geocentric embedding works better than using raw depth images for learning feature representations with convolutional neural networks. Our final object detection system achieves an average precision of 37.3%, which is a 56% relative improvement over existing methods. We then focus on the task of instance seg- mentation where we label pixels belonging to object instances found by our detector. For this task, we propose a decision forest approach that classifies pixels in the detection window as foreground or background using a family of unary and binary tests that query shape and geocentric pose features. Finally, we use the output from our object detectors in an existing superpixel classification framework for semantic scene segmentation and achieve a 24% relative improvement over current state-of-the-art for the object categories that we study. We believe advances such as those represented in this paper will facilitate the use of perception in fields like robotics.


#### Citing ####
If you find this code useful in your research, please consider citing:

    @incollection{guptaECCV14,
      author = {Saurabh Gupta and Ross Girshick and Pablo Arbelaez and Jitendra Malik},
      title = {Learning Rich Features from {RGB-D} Images for Object Detection and Segmentation},
      booktitle = ECCV,
      year = {2014},
    }

#### License ####
This code (RCNN-Depth) is released under the Simplified BSD License (refer to the LICENSE file for details). License files for individual packages are included within individual folders.

#### Installation Instructions ####
0. Create directory, checkout eccv14-code, utils, rgbdutils, nyu-hooks

  ```shell
  mkdir rcnn-depth && cd rcnn-depth
  git clone https://github.com/s-gupta/rcnn-depth.git eccv14-code
  git clone https://github.com/s-gupta/rgbdutils.git eccv14-code/rgbdutils
  git clone https://github.com/s-gupta/utils.git eccv14-code/utils
  git clone https://github.com/s-gupta/nyu-hooks.git eccv14-code/nyu-hooks
  ```

0. Checkout caffe-code 

  ```shell
  git clone https://github.com/BVLC/caffe.git eccv14-code/caffe
  cd eccv14-code/caffe
  git checkout e5cc609138a0bc4ce5177a67cf84952756d11b38
  cd ../../
    ```
  
0. Get the data (color image, depth images, rawdepth images, splits, ground truth, tasks), and external model data (Caffe trained Imagenet model, structured forests BSDS model). (Note: Some of these download links don't work with some versions of Google Chrome, please use wget or Firefox or Safari to download them.)

  ```shell
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-data.tgz
  tar -xf eccv14-data.tgz
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-external-data.tgz
  tar -xf eccv14-external-data.tgz
  ```

0. Get precomputed models.
 
  ```shell
  wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-models.tgz
  tar -xf eccv14-models.tgz 
  ```

0. [Optional] Get precomputed results

  ```shell
  wget ftp://ftp.cs.berkeley.edu/pub/projects/vision/rcnn-depth-eccv14/eccv14-cachedir.tgz
  tar -xf eccv14-cachedir.tgz # Also contains the models
  cd ..
  ```

#### Building ####
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
  cd ../../../
  ```
  
0. Build toolboxes, MCG, RCNN. Start MATLAB in the folder eccv14-code

  ```matlab
  mcg_build();
  rcnn_build();
  toolboxCompile();
  structured_edges_build();
  ```

#### Inference ####
##### Contour Detection, UCMs, Region Proposals and Detection on a new image #####
  ```matlab
  %% Should produce results on the image in demo-data 
  demo();
  
  %% See run_all.m, demo.m to figure out what units / data type
  % each image is in.
  [E, ucm2, candidates, detection_scores_no_nms, cls] = ...
    run_all(color_image, depth_image, rawdepth_image, camera_matrix, []);
  ```
  
#### Training ####
##### Contour Detection #####
0. Run the following in MATLAB (will require you to adapt jobParallel to run on your cluster / machine).

  ```matlab
  jobName = 'compute_edge_cues'; script_edges;
  jobName = 'train_edge_model'; script_edges;
  jobName = 'test_edge_model'; script_edges;
  ```

##### UCMs and Region Proposals #####
0. Run the following in MATLAB

  ```matlab
  jobName = 'edges_to_ucms'; script_regions;
  jobName = 'benchmark_multi_ucm'; script_regions;
  jobName = 'pareto'; script_regions;
  jobName = 'cache-mcg'; script_regions;
  jobName = 'rank_training'; script_regions;
  jobName = 'region-detect'; script_regions;
  ```

##### Object detectors: finetuning, training #####
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
  res = rcnn_all('task-detection', 'rgb_hha', 'trainval', 'test');
  ```
  
#### Notes ####
0. Although, the run_all function will produce outputs as produced by our system, if you are running the algorithm on a large set of images, it will be best to compute in batches using ``` scripts/script_regions.m``` and ```scripts/script_detection.m```
0. The current release of the code does not produce any instance segmentation or semantic segmentation. We will try to release code for those parts as soon as possible. If you are particularly interested in these 2 parts please let me know.
0. The current version of the code also does not include detectors finetuned on synthetic data, hence the performance is 1% worse than what we report in the paper (35.94% as oppsoed to 37.3%). Also, we accidently included 'cabinet' as a category in the finetuing, but never trained R-CNN detectors for the same. This version of the code mimicks the settings that were used to generate the tables in the papers.
0. For contour detection the provided model produces the solid blue line in Figure 3 in the paper.
0. We can also provide pre-computed results on the dataset. Please let me know if you want them.

#### FAQs ####
0. If you run out of memory on your GPU while extracting CNN features for object detection, you can decrease the batch size (which is currently set to 256). You will need to edit this in the following files: ```nyud2_finetuning / imagenet_*_256_fc6.prototxt``` and ```rcnn / rcnn_create_model.m```.

#### Acknowledgements ####
This work was sponsored by ONR SMARTS MURI N00014-09-1-1051, ONR MURI N00014-10-1-0933 and a Berkeley Fellowship. The GPUs used in this research were generously donated by the NVIDIA Corporation. We are also thankful to Bharath Hariharan, for all the useful discussions. We also thank Piotr Doll ́ar for helping us with their contour detection code.

#### Contact ####
If you find bugs / have questions, please let me know: Saurabh Gupta (sgupta [at] eecs [dot] berkeley [dot] edu)
