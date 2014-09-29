% ********************************************************************
% * mex File compiling code for Random Forest (for linux)
% * mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
% * Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
% * License: GPLv2
% * Version: 0.02 
% ********************************************************************/

function compile_linux

    % a simple way to compile on linux is to call the Makefile with 'make mex'
    system('rm *.mexglx *.mexa64;');

    %Matlab mex requires optimization to be all set in the mexopts.sh(or
    %.bat) file. So set it there not here
    
    %if you want to emulate what the makefile does ucomment below 2 lines:
    mex src/mex_regressionRF_train.cpp src/reg_RF.cpp src/cokus.cpp -o mexRF_train -DMATLAB 
    mex src/mex_regressionRF_predict.cpp src/reg_RF.cpp src/cokus.cpp -o mexRF_predict -DMATLAB 

%     system('make mex;');

    fprintf('Mex compiled\n')

