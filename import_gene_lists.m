function import_gene_lists(force)

if(~exist('force','var') || isempty(force))
    force = false;
end

% gene_info.txt from http://www.genenames.org/cgi-bin/download
fname1 = 'data/gene_info.txt';
fname2 = 'data/gene_info.mat';
if(~exist(fname2,'file') || force)
    websave(fname1,['http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=md_eg_id&col=md_ensembl_id&',...
        'status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit']);
    % get gene ids and symbols from txt file
    gene_info = readtable(fname1,'Delimiter','\t');
    gene_info.EnsemblGeneID = gene_info.EnsemblID_suppliedByEnsembl_;
    gene_info.EnsemblID_suppliedByEnsembl_ = [];
    gene_info.EntrezGeneID = gene_info.EntrezGeneID_suppliedByNCBI_;
    gene_info.EntrezGeneID_suppliedByNCBI_ = [];
    save(fname2,'gene_info'); 
end

% Housekeeping genes
fname1 = 'data/HK_genes.txt';
fname2 = 'data/HK_genes.mat';
if(~exist(fname2,'file') || force)
    websave(fname1, 'http://www.tau.ac.il/~elieis/HKG/HK_genes.txt');
    hk_genes = readtable(fname1,'Delimiter','space','ReadVariableNames',false);
    hk_genes = hk_genes.Var1;
    save(fname2,'hk_genes');
end

% LM22 gene set
fname1 = 'data/nmeth.3337-S2.xls';
fname2 = 'data/lm22.mat';
if(~exist(fname2,'file') || force)
    websave(fname1, 'http://www.nature.com/nmeth/journal/v12/n5/extref/nmeth.3337-S2.xls');
    lm22 = readtable(fname1,'Range','C14:C600');
    genes = unique(lm22.GeneSymbol);
    genes = genes(2:end);
    save(fname2,'genes');
end

% Table S3 gene set
fname1 = 'data/aad0501_Table_S3.xlsx';
fname2 = 'data/table_s3.mat';
if(~exist(fname2,'file') || force)
    websave(fname1, 'http://science.sciencemag.org/highwire/filestream/677381/field_highwire_adjunct_files/7/aad0501_Table_S3.xlsx');
    genes = readtable(fname1,'Range','A5:F200');
    genes = genes{:,:};
    genes = strrep(genes,'''','');
    genes = unique(genes(:));
    genes = genes(2:end);
    save(fname2,'genes');
end

% Table S12
fname1 = 'data/aad0501_Table_S12.xlsx';
fname2 = 'data/table_s12.mat';
if(~exist(fname2,'file') || force)
    websave(fname1, 'http://science.sciencemag.org/content/sci/suppl/2016/04/07/352.6282.189.DC1/aad0501_Table_S12.xlsx');
    genes = readtable(fname1,'Range','A7:A300');
    genes = genes{:,:};
    genes = unique(genes(:));
    genes = genes(2:end);
    save(fname2,'genes');
end
