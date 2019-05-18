function [Yhard, Yhard1,vv] = forestTest_new(model, X, opts)
    % X is NxD, where rows are data points
    % model comes from forestTrain()
    % Yhard are hard assignments to X's, Ysoft is NxK array of
    % probabilities, where there are K classes.
    
    if nargin<3, opts= struct; end
    numTrees= length(model.treeModels);
    u= model.treeModels{1}.classes; % Assume we have at least one tree!
    Ysoft= zeros(size(X,1), length(u));
    vv=[];
    for i=1:numTrees
        [~, ysoft] = treeTest(model.treeModels{i}, X, opts);
         Ysoft1 = ysoft/1;
         [~, ix]= max(Ysoft1, [], 2);
         Yhard(:,i) = u(ix);
        
         Ysoft= Ysoft + ysoft;
         vv=[vv Ysoft1(:,2)];
    end
    
    Ysoft = Ysoft/numTrees;
    [~, ix]= max(Ysoft, [], 2);
    Yhard1 = u(ix);
end
