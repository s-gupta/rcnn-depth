function rfImpute_test
    %Testing code, which uses the rfImpute_reg function below it

    %load the diabetes dataset
    load data/diabetes

    %modify so that training data is NxD and labels are Nx1, where N=#of
    %examples, D=# of features

    X = diabetes.x;
    Y = diabetes.y;

    [N D] =size(X);
    %randomly split into 250 examples for training and 50 for testing
    randvector = randperm(N);

    X_trn = X(randvector(1:250),:);
    Y_trn = Y(randvector(1:250));
    X_tst = X(randvector(251:end),:);
    Y_tst = Y(randvector(251:end));
    
    X_tst_orig = X_tst;
    X_tst(1,:)=NaN;
    indices_randomly_changed = find(isnan(X_tst));
    X_new = rfImpute_reg(X_tst,Y_tst);
    fprintf('Orig\tNew\n');
    [X_new(indices_randomly_changed)  X_tst_orig(indices_randomly_changed)]
    
    X = rand(100,2);X_orig = X;
    Y = 2*X(:,1)+3*X(:,2);
    X(1,:) = NaN;
    indices_randomly_changed = find(isnan(X));
    X_new = rfImpute_reg(X,Y);
    fprintf('Orig\tNew\n');
    [X_new(indices_randomly_changed)  X_orig(indices_randomly_changed)]
    
    
    
    X = [(1:50)' (51:100)'];
    Y = 2*X(:,1)+X(:,2);
    X(1:10,:) = NaN;
    indices_randomly_changed = find(isnan(X));
    X_new = rfImpute_reg(X,Y);
    fprintf('New\tOrig\n');
    [X_new(indices_randomly_changed)  X_orig(indices_randomly_changed)]
    1;
    
    
function X_out = rfImpute_reg(X_tst,Y_tst)
    iter = 50;
    ntree = 300;
    [X_out,hasNA] = na_roughfix(X_tst);
    
    
    for i=1:iter
        clear extra_options
        extra_options.proximity = 1; %(0 = (Default) Dont, 1=calculate)

        model = regRF_train(X_out,Y_tst, 100, 0, extra_options);
        Y_hat = regRF_predict(X_tst,model);
        %find proximity and then use it to fix values        
        prox = model.proximity;
        
        for j=1:length(hasNA)
            jj = hasNA(j);
            miss = find(isnan(X_tst(:,jj)));
            without_miss = 1:size(X_tst,1);
            without_miss(miss)=[];
            
            %have to add code for categorical variables
            %for non-categorical variables use the following
            
            sumprox = sum(prox(miss,without_miss)')';
            X_out(miss,jj)
            X_out(miss,jj) = prox(miss,without_miss) * X_out(without_miss,jj) ./  sumprox;
            X_out(miss,jj)
        end
    end
    


function [xout,hasNA] = na_roughfix(x)
    xout = x;
    hasNA=find(sum(isnan(x)));
    for i=1:length(hasNA)
        nan_indx = find(isnan(x(:,hasNA(i))));
        indx = 1:size(x,1);
        indx(nan_indx)=[];
        median(x(indx,hasNA(i)))
        xout(nan_indx,hasNA(i)) = median(x(indx,hasNA(i)));
    end
    find(sum(isnan(xout)))
    1;