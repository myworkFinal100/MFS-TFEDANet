function [model, weights] = trainMFSRegression(Xn, y, method)

% Normalize features
% mu = mean(X,1);
% sigma = std(X,[],1);
% sigma(sigma==0) = 1;
% Xn = (X - mu) ./ sigma;

switch upper(method)

    case 'RF'
        rng(1,'twister')
        model = fitrensemble(Xn, y, ...
            'Method','Bag', ...
            'NumLearningCycles',200, ...
            'Learners',templateTree('MinLeafSize',5));

        weights = predictorImportance(model);

    case 'LSBOOST'
        model = fitrensemble(Xn, y, ...
            'Method','LSBoost', ...
            'NumLearningCycles',300, ...
            'Learners',templateTree('MaxNumSplits',20));

        weights = predictorImportance(model);

    case 'GPR'
        model = fitrgp(Xn,y, ...
            'KernelFunction','ardsquaredexponential', ...
            'Standardize',false);

        weights = 1./model.KernelInformation.KernelParameters(1:end-1);

    case 'SVM'
        model = fitrsvm(Xn,y, ...
            'KernelFunction','gaussian', ...
            'Standardize',false);

        weights = ones(1,size(Xn,2));

    otherwise
        error('Unknown regression model.');

end

weights = weights ./ sum(weights);

end