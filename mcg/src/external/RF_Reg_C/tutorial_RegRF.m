% A simple tutorial file to interface with RF
% Options copied from http://cran.r-project.org/web/packages/randomForest/randomForest.pdf

%run plethora of tests
clear extra_options
clc
close all

%compile everything
if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
   compile_windows
else
   compile_linux
end

total_train_time=0;
total_test_time=0;

%diabetes
load data/diabetes
 
%modify so that training data is NxD and labels are Nx1, where N=#of
%examples, D=# of features

X = diabetes.x;
Y = diabetes.y;

[N D] =size(X);
%randomly split into 400 examples for training and 42 for testing
randvector = randperm(N);

X_trn = X(randvector(1:400),:);
Y_trn = Y(randvector(1:400));
X_tst = X(randvector(401:end),:);
Y_tst = Y(randvector(401:end));


 
% example 1:  simply use with the defaults
    model = regRF_train(X_trn,Y_trn);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 1: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
 
% example 2:  set to 100 trees
    model = regRF_train(X_trn,Y_trn, 100);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 2: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

% example 3:  set to 100 trees, mtry = 2
    model = regRF_train(X_trn,Y_trn, 100,2);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 3: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

% example 4:  set to defaults trees and mtry by specifying values as 0
    model = regRF_train(X_trn,Y_trn, 0, 0);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 4: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

% % example 5: set sampling without replacement (default is with replacement)
    extra_options.replace = 0 ;
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 5: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

% example 6: sampsize example
    %  extra_options.sampsize =  Size(s) of sample to draw. For classification, 
    %                   if sampsize is a vector of the length the number of strata, then sampling is stratified by strata, 
    %                   and the elements of sampsize indicate the numbers to be drawn from the strata.
    clear extra_options
    extra_options.sampsize = size(X_trn,1)*2/3;
    
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 6: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    
% example 7: nodesize
    %  extra_options.nodesize = Minimum size of terminal nodes. Setting this number larger causes smaller trees
    %                   to be grown (and thus take less time). Note that the default values are different
    %                   for classification (1) and regression (5).
    clear extra_options
    extra_options.nodesize = 7;
    
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 7: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
        

% example 8: calculating importance
    clear extra_options
    extra_options.importance = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 8: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    
    %model will have 3 variables for importance importanceSD and localImp
    %importance = a matrix with nclass + 2 (for classification) or two (for regression) columns.
    %           For classification, the first nclass columns are the class-specific measures
    %           computed as mean decrease in accuracy. The nclass + 1st column is the
    %           mean decrease in accuracy over all classes. The last column is the mean decrease
    %           in Gini index. For Regression, the first column is the mean decrease in
    %           accuracy and the second the mean decrease in MSE. If importance=FALSE,
    %           the last measure is still returned as a vector.
    figure('Name','Importance Plots')
    subplot(3,1,1);
    bar(model.importance(:,end-1));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Accuracy');
    
    subplot(3,1,2);
    bar(model.importance(:,end));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Gini index');
    
    
    %importanceSD = The ?standard errors? of the permutation-based importance measure. For classification,
    %           a D by nclass + 1 matrix corresponding to the first nclass + 1
    %           columns of the importance matrix. For regression, a length p vector.
    model.importanceSD
    subplot(3,1,3);
    bar(model.importanceSD);xlabel('feature');ylabel('magnitude');
    title('Std. errors of importance measure');

% example 9: calculating local importance
    %  extra_options.localImp = Should casewise importance measure be computed? (Setting this to TRUE will
    %                   override importance.)
    %localImp  = a D by N matrix containing the casewise importance measures, the [i,j] element
    %           of which is the importance of i-th variable on the j-th case. NULL if
    %          localImp=FALSE.
    clear extra_options
    extra_options.localImp = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 9: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

    model.localImp
    
% example 10: calculating proximity
    %  extra_options.proximity = Should proximity measure among the rows be calculated?
    clear extra_options
    extra_options.proximity = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 10: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

    model.proximity
    

% example 11: use only OOB for proximity
    %  extra_options.oob_prox = Should proximity be calculated only on 'out-of-bag' data?
    clear extra_options
    extra_options.proximity = 1; %(0 = (Default) Don't, 1=calculate)
    extra_options.oob_prox = 0; %(Default = 1 if proximity is enabled,  Don't 0)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 11: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));


% example 12: to see what is going on behind the scenes    
%  extra_options.do_trace = If set to TRUE, give a more verbose output as randomForest is run. If set to
%                   some integer, then running output is printed for every
%                   do_trace trees.
    clear extra_options
    extra_options.do_trace = 1; %(Default = 0)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 12: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

% example 13: to see what is going on behind the scenes    
%  extra_options.keep_inbag Should an n by ntree matrix be returned that keeps track of which samples are
%                   'in-bag' in which trees (but not how many times, if sampling with replacement)

    clear extra_options
    extra_options.keep_inbag = 1; %(Default = 0)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 13: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    
    model.inbag


% example 14: getting the OOB MSE rate. model will have mse field
    model = regRF_train(X_trn,Y_trn);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 14: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    
    figure('Name','OOB error rate');
    plot(model.mse); title('OOB MSE error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
    
% 
% example 15: nPerm
%               Number of times the OOB data are permuted per tree for assessing variable
%               importance. Number larger than 1 gives slightly more stable estimate, but not
%               very effective. Currently only implemented for regression.
    clear extra_options
    extra_options.importance=1;
    extra_options.nPerm = 1; %(Default = 0)
    model = regRF_train(X_trn,Y_trn,100,2,extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 15: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    
    figure('Name','Importance Plots nPerm=1')
    subplot(2,1,1);
    bar(model.importance(:,end-1));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Accuracy');
    
    subplot(2,1,2);
    bar(model.importance(:,end));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Gini index');
    
    %let's now run with nPerm=3
    clear extra_options
    extra_options.importance=1;
    extra_options.nPerm = 3; %(Default = 0)
    model = regRF_train(X_trn,Y_trn,100,2,extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 15: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    
    figure('Name','Importance Plots nPerm=3')
    subplot(2,1,1);
    bar(model.importance(:,end-1));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Accuracy');
    
    subplot(2,1,2);
    bar(model.importance(:,end));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Gini index');
    
% example 16: corr_bias (not recommended to use)
    clear extra_options
    extra_options.corr_bias=1;
    model = regRF_train(X_trn,Y_trn,100,2,extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 16: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    

% example 17: keep all tree predictions (returns prediction_per_tree a Ntst x Ntree matrix)
    clear extra_options
    model = regRF_train(X_trn,Y_trn,100,2);
    test_options.predict_all=1;
    [Y_hat,prediction_per_tree] = regRF_predict(X_tst,model,test_options);
    fprintf('\nexample 17: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));

% example 18: print after every tree. helpful when there is lots of data
% and tree creation is slow
    clear extra_options
    extra_options.print_verbose_tree_progression = 1; %(Default = 0)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 18: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));    

% example 19: get the nodes for test examples
%   terminal node indicators (ntest by ntree matrix)
    clear extra_options
    extra_options.nodes = 1;
    [Y_hat, tmp, nodes] = regRF_predict(X_tst,model,extra_options);
    fprintf('\nexample 19: node information now in node variable\n');

% example 20: using categories
    fprintf('\nexample 20: example of using categories\n');
    load data/diabetes
    X = diabetes.x;
    Y = diabetes.y;

    default_scale = 8; %X values will be scaled between 0-(default_scale-1)
    for i=1:size(X,2)
        min_0_X = double(X(:,i) - min(X(:,i))); % minimum value will be now 0.
        X_i_between_0_1 = min_0_X/max(min_0_X); %now the range will be 0-1
        X_i_between_0_scale_2_to_n = round(X_i_between_0_1 * (default_scale-1));
        X(:,i) = X_i_between_0_scale_2_to_n;
    end
    % X is fully categorized and ranges from values 0:7

    % construct the category information for X
    % .categories is a 1 x D vector saying what features are categorical and what are not categorical
    % true = categorical, false = numeric

    extra_options.categorical_feature = ones(1,size(X,2)); % all features are categorical

    % to choose features 1:5 as non-categorcal, use below code
    % extra_options.categorical_feature(1:5) = false; 

    % without categories
    model=regRF_train(X,Y, 500, 3);
    y_hat = regRF_predict(X, model);
    err_without_categories = sum(abs(y_hat-Y))
    avg_size_of_trees = mean(model.ndtree)

    % with categories
    model = regRF_train(X,Y, 500, 3, extra_options);
    y_hat = regRF_predict(X, model);
    err_with_categories = sum(abs(y_hat-Y))
    avg_size_of_trees = mean(model.ndtree)

    fprintf('\terror with/without categories should be similar, categorical trees might be smaller than numerical trees\n')
    