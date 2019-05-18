function [BestIndividual]=MultiobjectiveClustering(tsneX)
%load('data/tsne_data_8k_merged_Fast.mat')
M=2;
p1 = [19 13  7  5  4  3  3  2  3];
p2 = [ 0  0  0  0  1  2  2  2  2];
p1 = p1(M-1);
p2 = p2(M-1);
[N,W] = F_weight(p1,p2,M);
Generations=1;
MaxValue   = [1.5 20];
MinValue   = [0.5 10];
D=2;
Population = rand(N,D);
Population = Population.*repmat(MaxValue,N,1)+(1-Population).*repmat(MinValue,N,1)

W(W==0) = 0.000001;
T = 4;
%邻居判断
    B = zeros(N);
    for i = 1 : N
        for j = i : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
for i=1:1:N
Eps=Population(i,1);
MinPts=floor(Population(i,2));
cluster = DBSCAN(tsneX,Eps,MinPts);
cluster(cluster==-1) = 0; 
if sum(cluster)~=0
XX = matpy.mat2nparray(tsneX);
labels= matpy.mat2nparray(cluster);
FunctionValue(i,1)=-py.sklearn.metrics.silhouette_score(XX,labels, pyargs('metric','euclidean'));
FunctionValue(i,2)=-py.sklearn.metrics.calinski_harabaz_score(XX,labels);
else
FunctionValue(i,1)=1000;
FunctionValue(i,2)=1000;
end
end
Z = min(FunctionValue);

A = 1;
    Coding='Real';
    Boundary=[MaxValue;MinValue];
        %开始迭代
        for Gene = 1 : Generations
 for i = 1 : N
            %归一化
            Fmax = max(FunctionValue);
            Fmin = Z;
            FunctionValue = (FunctionValue-repmat(Fmin,N,1))./repmat(Fmax-Fmin,N,1);
            
            %选出父母
            k = randperm(T);
            k = B(i,k(1:2));

            %产生子代
            Offspring = P_generator([Population(k(1),:);Population(k(2),:)],Boundary,Coding,1);
           Eps=Offspring(1);
            MinPts=floor(Offspring(2));
            cluster = DBSCAN(tsneX,Eps,MinPts);
            cluster(cluster==-1) = 0; 
            if sum(cluster)~=0
                XX = matpy.mat2nparray(tsneX);
                labels= matpy.mat2nparray(cluster);
                OffFunValue(1)=-py.sklearn.metrics.silhouette_score(XX,labels, pyargs('metric','euclidean'));
                OffFunValue(2)=-py.sklearn.metrics.calinski_harabaz_score(XX,labels);
            else
                OffFunValue(1)=1000;
                OffFunValue(2)=1000;
            end
            OffFunValue = (OffFunValue-Fmin)./(Fmax-Fmin);
            
            %更新最优理想点
            Z = min(Z,OffFunValue);

            %更新邻居个体
            for j = 1 : T
                if A == 1
                    g_old = max(abs(FunctionValue(B(i,j),:)-Z).*W(B(i,j),:));
                    g_new = max(abs(OffFunValue-Z).*W(B(i,j),:));
                elseif A == 2
                    d1 = abs(sum((FunctionValue(B(i,j),:)-Z).*W(B(i,j),:)))/norm(W(B(i,j),:));
                    g_old = d1+5*norm(FunctionValue(B(i,j),:)-(Z+d1*W(B(i,j),:)/norm(W(B(i,j),:))));               
                    d1 = abs(sum((OffFunValue-Z).*W(B(i,j),:)))/norm(W(B(i,j),:));
                    g_new = d1+5*norm(OffFunValue-(Z+d1*W(B(i,j),:)/norm(W(B(i,j),:))));
                end
                if g_new < g_old
                    %更新当前向量的个体
                    Population(B(i,j),:) = Offspring;
                    FunctionValue(B(i,j),:) = OffFunValue;
                end
            end

            %反归一化
            FunctionValue = FunctionValue.*repmat(Fmax-Fmin,N,1)+repmat(Fmin,N,1);

 end

end
 FrontValue = P_sort(FunctionValue);
 BestIndividual=Population(FrontValue(1),:);