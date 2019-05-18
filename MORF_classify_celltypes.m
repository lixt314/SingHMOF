function [predicted_celltype,Testaccuracy] = MORF_classify_celltypes(mappedA,training_celltype);

rng(0);

%% train decision tree classifier 
tic;
id = training_celltype~=0;
id1=(training_celltype==0);
Index=find(training_celltype~=0)
Cell_1=find(training_celltype==1)

NP_1=randsample(Cell_1,floor(size(Cell_1,2)*0.1));
Cell_2=find(training_celltype==2)

NP_2=randsample(Cell_2,floor(size(Cell_2,2)*0.1));
Cell_3=find(training_celltype==3)

NP_3=randsample(Cell_3,floor(size(Cell_3,2)*0.1));
Cell_4=find(training_celltype==4)

NP_4=randsample(Cell_4,floor(size(Cell_4,2)*0.1));
Cell_5=find(training_celltype==5)

NP_5=randsample(Cell_5,floor(size(Cell_5,2)*0.1));
Cell_6=find(training_celltype==6)

NP_6=randsample(Cell_6,floor(size(Cell_6,2)*0.1));
Cell_7=find(training_celltype==7)

NP_7=randsample(Cell_7,floor(size(Cell_7,2)*0.1));
Cell_8=find(training_celltype==8)

NP_8=randsample(Cell_8,floor(size(Cell_8,2)*0.1));
Cell_9=find(training_celltype==9)

NP_9=randsample(Cell_9,floor(size(Cell_9,2)*0.5));
NP=[NP_1 NP_2 NP_3 NP_4 NP_5 NP_6 NP_7 NP_8 NP_9];

S_response = training_celltype(NP)';

S_predictors = mappedA(NP,:);


response = training_celltype(id)';

predictors = mappedA(id,:);


opts= struct;
opts.depth= 10;
opts.numTrees=200;
opts.numSplits= 200;
opts.verbose= true;
opts.classifierID= [2]; % weak learners to use. Can be an array for mix of weak learners too
PredicteClass=[];
 
mm= forestTrain(S_predictors, S_response,opts);

    
    
    D=200;
    M=2;
    p1 = [19 13  7  5  4  3  3  2  3];
    p2 = [ 0  0  0  0  1  2  2  2  2];
    p1 = p1(M-1);
    p2 = p2(M-1);
    [N,W] = F_weight(p1,p2,M);

    W(W==0) = 0.000001;
    T = 4;
    B = zeros(N);
    for i = 1 : N
        for j = i : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    MaxFES = 1000;
    Population = randi([0,1],N,D);

    for i=1:N
        if sum(Population(i,:)~=0)
    k=1;
        treeModels={};
        for j=1:1:D
            if Population(i,j)==1
                treeModels{k}= mm.treeModels{j};
                k=k+1;
            end
        end
      model.treeModels = treeModels;
      mx=model;
        
      [yhatTrain,Ysoft] = forestTest(mx, predictors);
        
      FunctionValue(i,1)=-size(find(yhatTrain==response),1)/size(response,1);
      FunctionValue(i,2)= sum(Population(i,:));
        else
            FunctionValue(i,1)=999999;
            FunctionValue(i,2)= 999999;    
        end
        % end
    end
    

    Z = min(FunctionValue)
    Coding='Binary';
    Boundary=[0;1];
    Generations=5
    for Gene = 1 : Generations
        [FrontValue,MaxFront] = P_sort(FunctionValue,'dd')
        mu=(N-FrontValue)./N;
        lambda = 1 - mu;
        for i = 1 : N
            Fmax = max(FunctionValue);
            Fmin = Z;         
            FunctionValue = (FunctionValue-repmat(Fmin,N,1))./repmat(Fmax-Fmin,N,1);
            for j = 1 : D
                if rand < lambda(i) % Should we immigrate?
                    % Yes - Pick a solution from which to emigrate (roulette wheel selection)
                    RandomNum = rand * sum(mu);
                    Select = mu(1);
                    SelectIndex = 1;
                    while (RandomNum > Select) && (SelectIndex < N)
                        SelectIndex = SelectIndex + 1;
                        Select = Select + mu(SelectIndex);
                    end
                    Offspring(j) = Population(SelectIndex, j); % this is the migration step
                else
                    Offspring(j) = Population(i, j); % no migration for this independent variable
                end

            end
            if sum(Offspring~=0)
                k=1;
                treeModels={};
            for j=1:1:D
                if Offspring(j)==1
                    treeModels{k}= mm.treeModels{j};
                    k=k+1;
                end
            end
                model.treeModels = treeModels;
                mx=model;    
                [yhatTrain,Ysoft] = forestTest(mx, predictors);
                OffFunValue(1)=-size(find(yhatTrain==response),1)/size(response,1);
                OffFunValue(2)= sum(Population(i,:));   
            else
                OffFunValue(1)=999999;
                OffFunValue(2)=999999;
            end
            OffFunValue = (OffFunValue-Fmin)./(Fmax-Fmin);
            Z = min(Z,OffFunValue);
            for j = 1 : T
                    g_old = max(abs(FunctionValue(B(i,j),:)-Z).*W(B(i,j),:));
                    g_new = max(abs(OffFunValue-Z).*W(B(i,j),:));
                if g_new < g_old
                    Population(B(i,j),:) = Offspring;
                    FunctionValue(B(i,j),:) = OffFunValue;
                end
            end
            FunctionValue = FunctionValue.*repmat(Fmax-Fmin,N,1)+repmat(Fmin,N,1);
        end       
        FunctionValue
    end

labellabel=yhatTrain';
New_labels=response';
size(labellabel);
size(New_labels);

for ijij=1:1:1
    Labels= labellabel(ijij,:);
    [position, Index1]=find(New_labels==1);
   
    TCell=size(find(New_labels(Index1)==labellabel(Index1)),2)./size(Index1,2)
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index3]=find(New_labels==2);
    BCell=size(find(New_labels(Index3)==labellabel(Index3)),2)./size(Index3,2) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index5]=find(New_labels==3);
    Macrophages=size(find(New_labels(Index5)==labellabel(Index5)),2)./size(Index5,2)   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index7]=find(New_labels==4);
    Dendritic=size(find(New_labels(Index7)==labellabel(Index7)),2)./size(Index7,2) 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index9]=find(New_labels==5);
    Naturalkillercells=size(find(New_labels(Index9)==labellabel(Index9)),2)./size(Index9,2)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index11]=find(New_labels==6);
    Endothelialcells=size(find(New_labels(Index11)==labellabel(Index11)),2)./size(Index11,2)     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index13]=find(New_labels==7);
    Cancer=size(find(New_labels(Index13)==labellabel(Index13)),2)./size(Index13,2)    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index15]=find(New_labels==8);
    OVarian=size(find(New_labels(Index15)==labellabel(Index15)),2)./size(Index15,2)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, Index17]=find(New_labels==9);
    Melanoma=size(find(New_labels(Index17)==labellabel(Index17)),2)./size(Index17,2)    
    
    Results=[TCell BCell Macrophages Dendritic Naturalkillercells Endothelialcells Cancer OVarian Melanoma]'

end
 for i=1:1:N
     k=1
     %treeModels=cell(1, sum(Population(i,:)));
     treeModels={};
     for j=1:1:D
        if Population(i,j)==1
            treeModels{k}= mm.treeModels{j};
            k=k+1;
        end
     end
      model.treeModels = treeModels;
      mx=model
      [predicted_celltype2,Ysoft] = forestTest(mx, mappedA);
      New_Ysoft=Ysoft(find(training_celltype==0),:);
      New_predicted_celltype2=predicted_celltype2(id1);
      
      New_predicted_celltype2(max(New_Ysoft,[],2)<0.99) = 0;
      predicted_celltype2(id1)=New_predicted_celltype2';
     %predicted_celltype2(id)=yhatTrain;
      predicted_celltype1(:,i)=predicted_celltype2';
      
    %Test_AUC(i)=roc_curve(Ysoft(:,2),test_y)
    
     Testaccuracy(i)=size(find(training_celltype'==predicted_celltype1(:,i)),1)/size(training_celltype',1)
  
 end
 [~,Index]=max(Testaccuracy)
 predicted_celltype=predicted_celltype1(:,Index);