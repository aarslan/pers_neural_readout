% classification() - run different classifcation algorithm for machine
%                   learning
%
% Usage: [YlatTr,YlatTest,B_new,best_lambda] = classification(Xtrain,Xtest,Ytrain,Ytest,classifier,optimization)
%
% Requierd input:     meaning (size or options) {default}
%
%   Xtrain          = matrix of data for training (Trials*features)
%   Xtest           = matrix of data for testing (Trials*features)
%   Ytrain          = colon vecotr of train labels (Trials)
%   Ytest           = colon vecotr of test labels (Trials)
%   classifier     = name of the algorithm to use ('fastLinearSVMAdrien'
%                     or 'L1')
%   optimization    = optimize the constraint parameter on the training set
%                    (on or off)
%
% Outputs:            meaning (size or options)
%
%   YlatTr          = colon vector of decimal value predicted by the
%                     algorithm for training set, to get the predicted value
%                     sign(YlatTr) (Trials)
%   YlatTest        = colon vector of decimal value predicted by the
%                     algorithm for test set, to get the predicted value
%                     sign(YlatTest) (Trials)
%   B_new           = weights (features)
%   best_lambda      = best constraint parameter get by optimization (x)
%
% Autors: Maxime, Serre Lab, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [YlatTr,YlatTest,B_new,best_lambda]= classification(Xtrain,Xtest,Ytrain,Ytest,classifier,optimization)

%==========================================================================
%                           Initializations
%==========================================================================

B_init = [];

% Add bias at the end
%-------------------
Xtrain      = [Xtrain ones(size(Xtrain,1),1)];
Xtest       = [Xtest ones(size(Xtest,1),1)];

%==========================================================================
%                           Classification
%==========================================================================

% Optimization
%----------------
if ~strcmp(optimization,'off')
    
    [best_lambda B_init] = parameterSearch(Xtrain,Ytrain,classifier);
else
    best_lambda = 1;
end

% Classification
%----------------
[YlatTr,YlatTest,B_new]= classif(Xtrain,Xtest,Ytrain,Ytest,classifier,best_lambda,B_init);

end





function [YlatTr,YlatTest,B_new]= classif(Xtrain,Xtest,Ytrain,Ytest,classifier,best_lambda,B_init)
% perform a classification of the input data using specified classifier
% classifiers need to be installed and in the path

switch classifier
    
    case 'fastLinearSVMAdrien'
        B_new = Lin_SVM_Keerthi(Xtrain,Ytrain,B_init,best_lambda);
        
    case 'L1'
        B_new = l1_ls(Xtrain,Ytrain,best_lambda);
        
    case 'logReg'
        out = logRegFun((Ytrain+1)/2,Xtrain',best_lambda);
        B_new = out.weights;
        
    case 'libSVM'
        model = svmtrain(Ytrain, Xtrain, ['-t 0 -c ', num2str(best_lambda)]);
        %B_new = -(model.SVs' * model.sv_coef);
        % to handle libSVM and liblinear "bug"
        % the B_new values can be sometimes reversed (see libSVM FAQ)
        labels = unique(Ytrain);
        if Ytrain(1) == labels(1);
            B_new = -(model.SVs' * model.sv_coef);
        elseif Ytrain(1) == labels(2);
            B_new = (model.SVs' * model.sv_coef);
        end
        
    case 'liblinear'
        model = train(Ytrain, sparse(double(Xtrain)), ['-s 2 -q -c ', num2str(best_lambda)]);
        %B_new = -model.w';
        % to handle libSVM and liblinear "bug" - the B_new values can be sometimes reversed (see libSVM FAQ)
        labels = unique(Ytrain);
        if Ytrain(1) == labels(1);
            B_new = -model.w';
        elseif Ytrain(1) == labels(2);
            B_new = model.w';
        end
end

% Value predicted
YlatTr              = Xtrain*B_new;
YlatTest            = Xtest*B_new;
B_new = B_new(1:end-1); % remove the bias to have clean weights
end



function [best_lambda B_init] = parameterSearch(Xtrain,Ytrain,classifier)
% ptraining = percentage of trial kept for training during the C optimization
% seb.crouzet@gmail.com

B_init      = [];
ptraining = 0.5;
ntraining = round((size(Xtrain,1)*ptraining)/2);
nfold = 8; % if 8 labs are open, there are thus done all at once (8 labs on the cluster)
foldPerf = nan(nfold,1);
switch classifier % specifies interval and grain for each classifier
    case 'libSVM'
        interval = [-15 4]; % C values | tested, seems to not be too bad
        grain = [1 0.1];
    case 'fastLinearSVMAdrien'
        interval = [-2 10]; 
        grain = [1 0.1];
    case 'logReg'
        interval = [3 5]; 
        grain = [0.2 0.02];
    case 'L1'
        interval = [-5 6]; 
        grain = [1 0.1];
    case 'liblinear'
        interval = [-15 4];
        grain = [1 0.1];
end

% grid search to find the best C
for pas = grain
    i=1;
    bestcv = 0;
    for log2c = interval(1):pas:interval(2)
        % repeat nfold times to be sure about the estimation for this C value
        % the nfolds can be parallelized if matlabpool is open, going 8 times faster if 8 labs are available
        parfor f=1:nfold
            train = [];
            labels = unique(Ytrain); % to handle both -1/1 and 0/1, maybe also multiclass
            for l = 1:length(labels)
                c = labels(l);
                j = find(Ytrain==c);
                t = randperm(length(j));
                train = cat(1, train, j(t(1:ntraining)));
            end
            test = setdiff(1:size(Ytrain,1), train)';
            [~,YlatTest,~] = classif(Xtrain(train,:),Xtrain(test,:),Ytrain(train),Ytrain(test),classifier,2^log2c,B_init);
            perf = measure_perf(YlatTest,Ytrain(test));
            foldPerf(f) = perf.accuracy;
        end
        fprintf('.')
        thisCVperf = mean(foldPerf);
        if (thisCVperf > bestcv),
            bestcv = thisCVperf;
            bestc = 2^log2c;
            bestlog2c = log2c;
        end
        i=i+1;
    end
    fprintf(1, '\n')
    interval(1) = bestlog2c-pas;
    interval(2) = bestlog2c+pas;
end
fprintf('\n')
best_lambda = bestc;
end