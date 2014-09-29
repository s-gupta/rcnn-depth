### Installation Instructions ###
0. Checkout eccv14-code, utils, rgbdutils, nyu-hooks
        ```sh
        git clone git@bitbucket.org:saurabhgupta/eccv14-code.git
        git clone git@bitbucket.org:saurabhgupta/rgbdutils.git eccv14-code/rgbdutils
        git clone git@bitbucket.org:saurabhgupta/utils.git eccv14-code/utils
        ```
0. Checkout caffe-code 
  <pre>
  $ git clone https://github.com/BVLC/caffe.git eccv14-code/caffe
  $ cd eccv14-code/caffe
  $ git checkout e5cc609138a0bc4ce5177a67cf84952756d11b38
  $ cd ../../
  </pre>
0. Get the data (color image, depth images, rawdepth images, splits, ground truth, tasks)
  <pre>
  $ wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-data.tgz
  $ tar -xf eccv14-data.tgz
  </pre>
0. Get the caffe models
  <pre>
  $ cd eccv14-code
  $ wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-caffe-data.tgz
  $ tar -xf eccv14-caffe-data.tgz 
  $ cd ..
  </pre>
0. Get precomputed models.
  <pre>
  $ cd eccv14-code
  $ wget http://www.cs.berkeley.edu/~sgupta/eccv14/eccv14-models.tgz
  $ tar -xf eccv14-models.tgz 
  $ cd ..
  </pre>

### Building ###
0. Build caffe (Adjust paths for CUDA / MATLAB in Makefile.config.example and copy to Makefile.config).
  <pre>
  >> cd eccv14-code/caffe
  >> make -j 16
  >> make -j 16 matcaffe
  </pre>
0. Build toolboxes, MCG, RCNN, ...
  <pre>
  >> mcg_build();
  >> rcnn_build();
  </pre>

### Running ###
#### Detection on a new image ####
  <pre>
  >> 
  </pre>
#### Edges, UCMs and Region Proposals ####
        ```matlab
        >> jobName = 'edges_to_ucms'; script_regions;
        >> jobName = 'benchmark_multi_ucm'; script_regions;
        >> jobName = 'pareto'; script_regions;
        >> jobName = 'cache-mcg'; script_regions;
        >> jobName = 'rank_training'; script_regions;
        >> jobName = 'region-detect'; script_regions;
        ```

#### Detector finetuning, training, testing ####
0. Run the following on the MATLAB prompt. At the end of this, you will get 2 finetuning commands that you need to use with caffe for finetuning.
  <pre>
  >> jobName = 'save_color'; script_detection;
  >> jobName = 'save_hha'; script_detection;
  >> jobName = 'write_window_file'; script_detection;
  </pre>

0. Use the finetuning commands that print out to finetune the CNN. Use the following to extract features and train the RCNN model.
  <pre>
  >> jobName = 'hha_cache_features'; script_detection;
  >> jobName = 'color_cache_features'; script_detection;
  >> res = rcnn_all('task-guptaetal', 'train', 'val');
  </pre>
