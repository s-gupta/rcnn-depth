%**************************************************************
%* mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
%* Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
%* License: GPLv2
%* Version: 0.02
%
% Calls Regression Random Forest
% A wrapper matlab file that calls the mex file
% This does prediction given the data and the model file
%**************************************************************

function [Y_hat,prediction_per_tree, nodes] = regRF_predict(X,model,extra_options)
    %function Y_hat = regRF_predict(X,model)
    %requires 2 arguments
    %X: data matrix
    %model: generated via regRF_train function
	if nargin<2
		error('need atleast 2 parameters,X matrix and model');
    end
    if exist('extra_options','var')
        if isfield(extra_options,'predict_all') 
            predict_all = extra_options.predict_all;
        end
        if isfield(extra_options,'nodes') 
            nodes = extra_options.nodes;
        end
    end
    
    if ~exist('predict_all','var'); predict_all=0;end
    
    if ~exist('nodes','var'); nodes=0;end
    
    if isfield(model, 'categorical_feature')
        % have to map prediction array to the correct categories
        for i=1:size(X,2)
            if model.categorical_feature(i) 
                tmp_uniques_in_feature = model.orig_uniques_in_feature{i};
                tmp_mapped_uniques_in_feature = model.mapped_uniques_in_feature{i};
                X_loc = X(:,i); %cannot change the original array which may cause chained change of categories to something totally wrong
                for j=1:length(tmp_uniques_in_feature)
                    indices_to_change = find( X(:,i) == tmp_uniques_in_feature(j) );
                    X_loc(indices_to_change) = tmp_mapped_uniques_in_feature(j);
                end
                X(:,i) = X_loc;
            end
        end
        ncat = model.ncat;
    else
        ncat = ones(1,size(X,2));
    end
    
    maxcat = max(ncat);
	[Y_hat,prediction_per_tree, nodes] = mexRF_predict(X',model.lDau,model.rDau,model.nodestatus,model.nrnodes,model.upper,model.avnode,model.mbest,model.ndtree,model.ntree,predict_all, nodes, int32(ncat), maxcat );
    
    if ~isempty(find(model.coef)) %for bias corr
        Y_hat = model.coef(1) + model.coef(2)*Y_hat;
    end
	clear mexRF_predict
    