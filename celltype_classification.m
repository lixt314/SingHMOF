function celltype_classification(pbmc_dataset,geneset,savefigure)
addpath('tools')
addpath('tools/bhtsne/')
close all


%%
fname = sprintf('data/normalized_data_%s.mat', pbmc_dataset);
if(~exist(fname,'file'))
    %% load pbmc data
    load(sprintf('data/Liver_%s.mat', pbmc_dataset));
    pbmc_data = data; %#ok<*NODEF>
    pbmc_sample = data.Properties.VariableNames;
    pbmc_sample = strrep(pbmc_sample,'P1TLH_[ATGC]*_','P1TLH');
    pbmc_sample = strrep(pbmc_sample,'P5TLH_','');
    pbmc_sample = strrep(pbmc_sample,'1','1_4k');
    
    tmp = regexp(pbmc_sample,'(?<=[ATGC]*_[0-9]_).*','match');
    pbmc_sample = [tmp{:}];
 
 
    pbmc_sample = strcat('pbmc_',pbmc_sample);
    
    clear data;
    %% merge and normalize
    [pbmc_data] =...
        normalizeData(pbmc_data);
    data = [pbmc_data];
    sample_id = [pbmc_sample];
    save(fname,'data','sample_id','-v7.3');
else
    load(fname);
end

fname = sprintf('data/classified_data_%s_%s.mat', pbmc_dataset, geneset);
if (~exist(fname,'file'))
    %% restrict to 'merged' gene set
    Xn=table2array(data);
    genelist=data.Properties.RowNames   
    [T]=sc_hvg(Xn,genelist,true,true);
    % Highly variable genes (HVGenes), FDR<0.05
    HVGenes=T.genes(T.fdr<0.05)
    fname12 = sprintf('data/Inf_genes.mat');
    if(~exist(fname12,'file'))
     [~,ib] = intersect(genelist, HVGenes);
     X=data{ib,:}';
     [RANKED, WEIGHT] = infFS( X,[],0.5,0,0);
     Inf_Index=ib(RANKED(1:300));
     Inf_genes=data.Properties.RowNames(Inf_Index);
     save(fname12,'Inf_genes');
    else
        load(fname12)
    end

    genes = data.Properties.RowNames;
    
    markers = load('data/marker_genes1.mat');
        
    selected_genes=markers.genes;

    selected_genes = unique(union(selected_genes, Inf_genes(1:100)));

    [~,ia] = intersect(genes, selected_genes); 
    save('selected_genes.mat','selected_genes')



    fname2 = sprintf('data/tsne_data_%s_%s_Fast.mat', pbmc_dataset, geneset);
    if(~exist(fname2,'file'))
         addpath('/tools/FIt-SNE-master');
         tsneX = fast_tsne(data{ia,:}');
        save(fname2,'tsneX');
    else
        load(fname2);
    end
   
    %%%%%%%%%%%%%%%%%%The paramemters obtained from evolutionary
    %%%%%%%%%%%%%%%%%%multiobjective clustering
%     epsilon=1.5;
%     ncells=18;
%     cluster = DBSCAN(tsneX,epsilon,ncells);
%     
%     cluster(cluster==-1) = 0;

    [BestCluster]=MultiobjectiveClustering(tsneX)
    epsilon=BestCluster(1);
    ncells=BestCluster(2);
     cluster = DBSCAN(tsneX,epsilon,ncells);
%     
     cluster(cluster==-1) = 0;

    load('colors.mat');


    dbscan_cols = parula(length(unique(cluster)));
    rng(0);
    dbscan_cols = dbscan_cols(randperm(max(unique(cluster))),:);
    dbscan_cols = [0.8*ones(1,3); dbscan_cols];

    fig = figure(1);
    clf
    lims = [min(tsneX(:)) max(tsneX(:))];
    for ci=1:length(indication)
        id = ~cellfun(@isempty,strfind(sample_id,indication{ci}));
        my_gscatter(tsneX(id,1),tsneX(id,2),cluster(id),dbscan_cols,symbs(ci),false,[],25,lims);
    end
    my_gscatter(tsneX(:,1),tsneX(:,2),cluster,[],{'none'},true,[],25,lims);

    fig.Color = 'w';
    fig.Position = [0 0 900 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);

    fname = sprintf('figures/dbscan_clusters_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
    end

    %% marker genes
    fig = figure(2);
    clf
    [celltype_matrix, celltype_expression, cellnames] = ...
        celltype_markers(data,tsneX,[2:10],0.3); 
    fig.Color = 'w';
    fig.Position = [0 0 800 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
    fname = sprintf('figures/marker_expression_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
    end

    %% training data
    training_celltype = assign_clusters(cluster,celltype_matrix,celltype_expression);


    fprintf(1,'%d/%d celltypes found for classification training\n',sum(unique(training_celltype)~=0),9);


    %% plot training set
    fig = figure(3);
    clf;

    lims = [min(tsneX(:)) max(tsneX(:))];
    for ci=1:length(indication)
        id = ~cellfun(@isempty,strfind(sample_id,indication{ci}));
        my_gscatter(tsneX(id,1),tsneX(id,2),training_celltype(id),cols,symbs(ci),false,[],25,lims);
    end
   legend('Unknown','hepatocytes cells','T cells','Mature B cells','Plasma cells','Endothelial cells','Inflammatory monocyte/macrophages cells',...
            'Non-inflammatory macrophages cells','Cholangiocyte cells','Stellate cells')

    fig.Color = 'w';
    fig.Position = [0 0 900 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
    fname = sprintf('figures/celltype_training_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
    end

    %% classify
    % Train deep autoencoder to map data
   fname5 = sprintf('data/Deep_data_%s_%s.mat', pbmc_dataset, geneset);
    if(~exist(fname5,'file'))
        no_dims=10;
        A=data{ia,:};
        A=A';
        addpath(genpath('E:\XiangtaoPaper\Bioinformatics\drtool'));
        layers = [ceil(size(A, 2) * 1.2) + 5 max([ceil(size(A, 2) / 4) no_dims + 2]) + 3 max([ceil(size(A, 2) / 10) no_dims + 1]) no_dims]
        disp(['Network size: ' num2str(layers)]);
        [network, mappedA] = train_deep_autoenc(A, layers, 0);
        save(fname5,'mappedA');
    else
        load(fname5);
    end

    [predicted_celltype,Testaccuracy] = MORF_classify_celltypes(mappedA,training_celltype);
    predicted_celltype=predicted_celltype';
    Testaccuracy=size(find(training_celltype==predicted_celltype),2)/size(training_celltype,2)
    %% plot classification result
    fig = figure(4);
    clf;
    lims = [min(tsneX(:)) max(tsneX(:))];
    for ci=1:length(indication)
        id = ~cellfun(@isempty,strfind(sample_id,indication{ci}));
        my_gscatter(tsneX(id,1),tsneX(id,2),predicted_celltype(id),cols,symbs(ci),false,[],25,lims);
    end

    fig.Color = 'w';
    fig.Position = [0 0 900 800];
    fig.PaperPositionMode = 'auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.PaperPosition(3:4);
    fname = sprintf('figures/celltype_classified_%s_%s', pbmc_dataset, geneset);
    if(savefigure)
        print(fig,'-dpng',[fname '.png'],'-r300');
        saveas(fig,[fname '.fig']);
    end

    fname3 = sprintf('data/classified_data_%s_%s.mat', pbmc_dataset, geneset);
    save(fname3,...
        'tsneX','cluster',...
        'sample_id','predicted_celltype','cellnames',...
        '-v7.3');

end

