function [data1] = normalizeData(data1)

    %if(~istable(data1))
    %    error('input data must be table')
   % end

    % load house-keeping genes
   if(~istable(data1))
        error('input data must be table')
    end
     
    % load house-keeping genes
    load('data/HK_genes.mat');

    % find common genes
    [~, ~, ia] = intersect(hk_genes,data1.Properties.RowNames);
  
    % convert to TPM scale
    data1 = logTrafo(data1,-1);

    % normalize to house-keeping gene expression
    hk_expr = mean(data1{ia,:},1);

    data1{:,:} = bsxfun(@times,data1{:,:},mean(hk_expr)./hk_expr);

    % convert back to log scale
    data1 = logTrafo(data1,1);
    %data3=data1;
    % find common genes
    %[~, ~, ia] = intersect(hk_genes,data1.Properties.RowNames);

%     % convert to TPM scale
%     %data1 = logTrafo(data1,-1);
 
%      [Xs]=sc_norm(A,'type','libsize');
%      Xs=log(Xs+1);
       
%      %Xs=run_magic(Xs);
% %     %figure; imagesc(Xs); title('Magic Normalized'); colorbar;
%      data2 = array2table(Xs);
%      data2.Properties.RowNames = data1.Properties.RowNames;
%      data2.Properties.VariableNames = data1.Properties.VariableNames;
%      data1=data2;
    
    % normalize to house-keeping gene expression
    %hk_expr = mean(data1{ia,:},1);

    %data1{:,:} = bsxfun(@times,data1{:,:},mean(hk_expr)./hk_expr);

    % convert back to log scale
    

end
