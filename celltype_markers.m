function [celltype_matrix, celltype_expression, cellnames, celltype] = celltype_markers(data,tsneX,id_celltype,thresh,id)

genes = {'Unknown',{},{},{};...
    'hepatocytes cells',{'HPR','GSTA2HPR','GSTA2','AKR1C1','SCD', 'HMGCS1','ACSS2',...
    'TM7SF2', 'SEC16B', 'SLBP','RND3', 'ARG1','BCHE', 'G6PC', 'GHR','ALDH6A1',...
    'RCAN1', 'AR','RPP25L', 'HSD11B1','HAMP', 'GHR', 'APOM'},{},{};...
    'T cells',{'CD2', 'TRAC', 'GZMK','CD3D', 'CD247', 'TRDC','STMN1', 'MKI67', ...
    'NKG2A','CD7', 'KLRB1', 'KLRD1', 'NKG7','CD69', 'NCAM1'},{},{};...
    'Mature B cells',{'CD20', 'LTB', 'CD37', 'CD79B'},{},{};...
    'Plasma cells',{'IGLC2', 'IGHG1', 'IGKC', 'CD79A', 'CD38'},{},{};...
    'Endothelial cells',{'MGP', 'SPARCL1', 'TMSF1', 'CLEC14A','IDI', 'IGFBP7','CTGF', 'VWF','CD32B', ...
    'LYVE1','STAB2','CCL14','CLEC1B','FCN2','S100A13', 'CD105', 'PECAM','PRSS23', 'RAMP3','INMT'},{},{};...
    'Inflammatory monocyte/macrophages cells',{'CD68', 'S100A8/9', 'LYZ', 'VCAN','IL18','CD74'},{},{};...
    'Non-inflammatory macrophages cells',{'CD68', 'CD5L','MARCO', 'VSIG4','CPVL','CD163', 'VCAM1'},{},{};...
    'Cholangiocyte cells',{'KRT19','EPCAM','FXDY2','CLDN4','CLDN10','MMP7','CXCL1','KRT7'},{},{};...
    'Stellate cells',{'ACTA2', 'COL1A1', 'COL1A2','TAGLN', 'COL3A1','DCN'},{},{};...
    'Sinusoidal endothelial cells (periportal)',{'MGP', 'SPARCL1', 'TMSF1', 'CLEC14A','IDI', 'IGFBP7','CTGF', 'VWF'},{},{};...
    'Sinusoidal endothelial cells (Central venous by IF)',{'CD32B','LYVE1','STAB2','CCL14','CLEC1B','FCN2','S100A13'},{},{};...
    'Portal endothelial cell',{'CD105', 'PECAM','PRSS23', 'RAMP3','INMT'},{},{}}; % J Clin Aesthet Dermatol. 2014 Jun; 7(6): 13?24.

if(~exist('id_celltype','var') || isempty(id_celltype))
    id_celltype = 1:size(genes,1);
end
if(~exist('thresh','var') || isempty(thresh))
    thresh = 0;
end
if(~exist('id','var') || isempty(id))
    id = true(1,size(tsneX,1));
end

lb = min(tsneX);
ub = max(tsneX);

%% reduce to id
data = data(:,id);
tsneX = tsneX(id,:);

%%
celltype_matrix = zeros(size(genes,1),size(data,2));
celltype_expression = zeros(size(genes,1),size(data,2));
celltype = zeros(1,size(data,2));
cellnames = genes(:,1)';

if length(id_celltype)<10
    m = 3;
else
    m = 4;
end

for i=1:length(id_celltype)
    j = id_celltype(i);
    
    [~, idx] = intersect(data.Properties.RowNames,genes{j,2});
    or_marker_expr = data{idx,:};
    or_marker_expr = bsxfun(@rdivide,...
        or_marker_expr,...
        max(or_marker_expr,[],2)); 
    or_marker_expr = mean(or_marker_expr,1);
    
    if ~isempty(genes{j,3})
        [~, idx] = intersect(data.Properties.RowNames,genes{j,3});
        and_marker_expr = data{idx,:};
        and_marker_expr = bsxfun(@rdivide,...
            and_marker_expr,...
            max(and_marker_expr,[],2));
        and_marker_expr = mean(and_marker_expr,1);
    else
        and_marker_expr = ones(size(or_marker_expr));
    end
    
    if ~isempty(genes{j,4})
        [~, idx] = intersect(data.Properties.RowNames,genes{j,4});
        not_marker_expr = data{idx,:};
        not_marker_expr = bsxfun(@rdivide,...
            not_marker_expr,...
            max(not_marker_expr,[],2));

        not_marker_expr = 1-mean(not_marker_expr,1);
    else
        not_marker_expr = ones(size(or_marker_expr));
    end
    
    marker_expr = or_marker_expr.*and_marker_expr.*not_marker_expr;
    
    marker_expr = bsxfun(@rdivide,marker_expr,max(marker_expr,[],2));
    
    celltype_matrix(j,:) = marker_expr>thresh;
    celltype_expression(j,:) = marker_expr;
    celltype(marker_expr>thresh) = j-1;
    
    [marker_expr, id_sort] = sort(marker_expr);
    
    subaxis(m,m,i,'S',0.02,'SV',0.03,'M',0.03,'MT',0.05)
    s = scatter(tsneX(id_sort,1),tsneX(id_sort,2),5,marker_expr);
    s.MarkerFaceColor = s.MarkerEdgeColor;
    s.MarkerFaceAlpha = 0.5;
    s.MarkerEdgeAlpha = 0.5;
    ax = gca;
    colormap grayorangered;
    ax.XTick = [];
    ax.YTick = [];
    xlim([-1.05*abs(lb(1)) 1.05*abs(ub(1))]);
    ylim([-1.05*abs(lb(2)) 1.05*abs(ub(2))]);
    box on;
    title(genes{j,1},'Fontname', 'Consolas')    
end
