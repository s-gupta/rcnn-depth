% ------------------------------------------------------------------------ 
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
% 
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2013
% ------------------------------------------------------------------------ 

FIRST INSTALL
- Change root_dir to point to the MCG root folder (the one with the file root_dir.m, instal.m, etc.)
- Change datasets/database_root_dir to point to your PASCAL2012 folder (the one with subfolders ImageSets, JPEGImages, etc.)
- Run install.m from the root dir to add the needed paths and do some checks
- If you need to re-build the library (the script install.m will tell if needed), run build.m

USAGE INSTALL
- Each time you restart your matlab, run install.m
- If you want to avoid this, add the paths permanently

GENERALITIES
- Three train subdivisions are provided, to validate the algorithm, if you want other subdivisions, run src/aux/create_train_samples.m
- Some partial results are stored automatically when computed. If loaded, the code will generally say it. If you need it to be recomputed, erase or rename that file.

UCMs AT MULTIPLE SCALES
- The technique is based on computing UCMs at different scales. The first step is therefore to compute them using "scripts_training/compute_all_ucms"
- If you prefer, you can download the pre-computed UCMs at multiple scales from the MCG website and put them inside "datasets/pascal2012", in a folder named "sf_mUCM"

PARETO LEARNING (See "scripts_trining/pareto_choose_point.m")
1) It runs pareto_learning.m (it takes a few minutes and stores the results),
2) it allows to choose the working point,
3) and it stores the selected parameters to file.

RANK TRAINING (See "scripts_trining/rank_training.m")
1) It gathers a random sample of examples from the selected pareto front point, computes their features and stores it (it takes a few hours).
2) It trains the ranking regressor. (it takes a few minutes)
- If you just want to re-train the regressor, the features will be loaded from file.
- If you want to re-sample the features, erase the stored file.

MCG CANDIDATES (See "scripts_trining/compute_mcg_cands.m")
- It computes the MCG candidates and their bounding boxes and stores them at files.
- It is prepared to run in parallel in a cluster: just run as many processes as wanted. Each process stores an empty file when it starts an image in particular. When it ends, it substitutes the file by the actual result.

EVALUATION
- See "scripts/benchmark_results.m" and ucomment the line that evaluates the current result 
