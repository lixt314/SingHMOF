function import_Liver_data(force)
%,data,matrix,genes,barcodes
addpath('tools/')

if(~exist('force','var') || isempty(force))
    force = false;
end

if(~exist('data/Liver_8k.mat','file') || force)
    
    data=readtable('/public/home/zxl/Research/Estimate/mouse/GSE115469_Data.csv','ReadRowNames',true);
    data = logTrafo(data, 1);
    save('data/Liver_8k.mat','data','-v7.3');
    clear;
    %data1 = readtable('/public/home/zxl/Research/Estimate/mouse/GSE115469_Data.csv','Delimiter','\t')
%     matrix = data1.data;
%     genes = readtable('data/genes.xls',...
%         'Filetype','Text','Delimiter','\t','ReadVariableNames',false)
%  
%     barcodes = readtable('data/barcodes.tsv',...
%         'Filetype','Text','Delimiter','\t','ReadVariableNames',false);
%     data = array2table(matrix,...
%         'RowNames',genes,...
%         'VariableNames',strrep(barcodes,'-','_'));
%     clear matrix;
end