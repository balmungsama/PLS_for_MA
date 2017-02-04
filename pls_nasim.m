
function [avg_ZSalience_X,avg_ZSalience_Y,pred_scores_X, pred_scores_Y,pls_out] = pls_loo(X,Y,var_norm);
% .
% .
% =========================================================================
% pls_loo:   Partial Least Squares (PLS) in a leave-one-out
% framework 
% 
% =========================================================================
%
% Syntax:
%              [avg_ZSalience_X,avg_ZSalience_Y,pred_scores_X, pred_scores_Y,pls_out]   =  pls_loo(X,Y,var_norm);
% Input: 
%         X,Y        = 2D data matrices, being compared. Should be formatted as:
%                          X = [N samples x Px variables] matrix of neuroimaging data
%                          Y = [N samples x Py variables] matrix of behavioural data
%                      
%
%        
%         var_norm       = mean centering/normalization method applied to both X
%                           and Y matrices  
%   
%                        0= no normalization, directly use input values 
%                        1= (column wise) mean centring  of X and Y 
%                        2= zscore X and Y 
%                        
%
% -------------------------------------------------------------------------
% Output:   "pls_out" list array of structures,includes results for every iteration of the leave-one-out  
%         
%             avg_ZSalience_X,avg_ZSalience_Y = bootstrap ratio for PLS components, averaged across leave-one-out iterations  
%             
%             
%            pred_scores_X, pred_scores_Y  = predicted subject scores, calculated by projecting the left out subject on the 
%                                             PLS component 
%            
%
%         pls_out.Salience_X, pls_out.Salience_Y    =  PLS components for X and Y,
%                                       respectively 
% 
%         pls_out.ZSalience_X, pls_out.ZSalience_Y  =  bootstrap ratio for PLS components for X and Y,
%                                       respectively 
%         pls_out.latent_X, pls_out.latent_Y        =  projection of X and Y data onto the PLS components 
%                                        (i.e scores) for N-1 subjects/samples from
%                                        which the PLS componsnts are
%                                        calculated 
%
%         pls_out.latent_X, pls_out.latent_Y        = projection of the left-out
%                                     subject/sample onto the PLS components
%   
%                        
%-------------------------------------------------------------------------
% 
% Nasim Shams
% version 13/01/2016
%
% edit 06/06/2016
% in PLS_between_MEG_fMRI function
% lines related to region analyais removed
% M,F,M_T, F_T outputs removed .
%%
nsub = size(Y,1);
ncomp = min(size(Y,2),size(X,2));
Xfull = X;
Yfull = Y;


for ii = 1 : nsub
    x = Xfull;
    y = Yfull;
    xo = Xfull(ii,:);
    yo = Yfull(ii,:);
    XO(ii,:)=xo;
    YO(ii,:)=yo;
    x(ii,:)=[];
    y(ii,:)=[];
    xm = mean(x);
    xstd = std(x);
    ym = mean(y);
    ystd= std(y);
    
    switch var_norm
        
        case 1
            x = detrend(x,'constant');
            y = detrend(y,'constant');
            xo = xo-xm;
            yo = yo-ym;
        case 2
            
            xo = (xo-xm)./xstd;
            yo = (yo-ym)./ystd;
            x = zscore(x);
            y = zscore(y);
            
    end   


    [pls_res(ii).Salience_X,pls_res(ii).Salience_Y,pls_res(ii).latent_X,pls_res(ii).latent_Y,pls_res(ii).ZSalience_X,pls_res(ii).ZSalience_Y,pls_res(ii).Sig_prob,pls_res(ii).Dcorr_pct] = PLS_between_MEG_fMRI(x,y,1,1);

    pls_res(ii).latent_Yo  = yo*pls_res(ii).Salience_Y;
    pls_res(ii).latent_Xo  = xo*pls_res(ii).Salience_X;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end


switch var_norm
    
    case 1
        Y = detrend(Y,'constant');
        X = detrend(X,'constant');
    case 2
        X = zscore(X);
        Y = zscore(Y);
        
end
  
% % % 

[pls_full.Salience_X,pls_full.Salience_Y,pls_full.latent_X,pls_full.latent_Y,pls_full.ZSalience_X,pls_full.ZSalience_Y,pls_full.Sig_prob,pls_full.Dcorr_pct] = PLS_between_MEG_fMRI(X,Y,1,1);

pls_sort = pls_res;
for ii = 1: nsub 
    [ind(ii,:),sg(ii,:)] = sort_eigen_images(pls_full.Salience_X,pls_res(ii).Salience_X);
    pls_sort(ii).ZSalience_X = bsxfun(@times,pls_res(ii).ZSalience_X(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).Salience_Y = bsxfun(@times,pls_res(ii).Salience_Y(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).Salience_X = bsxfun(@times,pls_res(ii).Salience_X(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).ZSalience_Y = bsxfun(@times,pls_res(ii).ZSalience_Y(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).latent_Y =  bsxfun(@times,pls_res(ii).latent_Y(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).latent_X =  bsxfun(@times,pls_res(ii).latent_X(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).latent_Yo =  bsxfun(@times,pls_res(ii).latent_Yo(:,ind(ii,:)),sg(ii,:));
    pls_sort(ii).latent_Xo =  bsxfun(@times,pls_res(ii).latent_Xo(:,ind(ii,:)),sg(ii,:));
   
end

pls_out = pls_sort;

for ii = 1: nsub 
    avg_ZSalience_X(:,:,ii) = pls_sort.ZSalience_X;
    avg_ZSalience_Y(:,:,ii) = pls_sort.ZSalience_Y;
end
avg_ZSalience_X = mean(avg_ZSalience_X,3);
avg_ZSalience_Y = mean(avg_ZSalience_Y,3);

pred_scores_X = reshape([pls_sort.latent_Xo],[ncomp,nsub]);
pred_scores_Y = reshape([pls_sort.latent_Yo],[ncomp,nsub]);
pred_corvals = diag(corr(pred_scores_X',pred_scores_Y'));

function [Salience_M,Salience_f,latent_M,latent_f,ZSalience_M,ZSalience_f,Sig_prob,Dcorr_pct,VSalience_M,VSalience_f] = PLS_between_MEG_fMRI(MEG_Map,BOLD_Map,num_regions,pvalue)

[Um,Dm,Vm] = svd(MEG_Map,'econ');
[Uf,Df,Vf] = svd(BOLD_Map,'econ');

[Umixed,Dmixed,Vmixed] = svd(Dm'*Um'*Uf*Df,'econ');
Salience_M =Vm*Umixed;
Salience_f = Vf*Vmixed;
clear Vm Vf
Dcorr = diag(Dmixed);
Dcorr_pct = round((Dcorr./sum(Dcorr)).*100);

for ires = 1:1000
    iperm = randperm(size(MEG_Map,1)); 
    [Umixed,Dmixed,Vmixed] = svd(Dm'*Um(iperm,:)'*Uf*Df,'econ');

    Dcorrp(:,ires) = diag(Dmixed);
end
Sig_prob = mean(bsxfun(@gt,Dcorr,Dcorrp),2);
Sig_Index = Sig_prob>=(1-pvalue);
% sum(Sig_Index)
% Sig_prob

Salience_M = Salience_M(:,Sig_Index);
Salience_f = Salience_f(:,Sig_Index);

latent_M   = MEG_Map*Salience_M;
latent_f   = BOLD_Map*Salience_f;
RSalience_M = 0;
RSalience_f = 0;
MSalience_M = 0;
MSalience_f = 0;

vRSalience_M = 0;
vRSalience_f = 0;
vMSalience_M = 0;
vMSalience_f = 0;


for bs = 1:1000
    %bs
    
    isub = ceil(size(MEG_Map,1)*rand(1,size(MEG_Map,1)));
    
    [bUm,bDm,bVm] = svd(MEG_Map(isub,:),'econ');
    [bUf,bDf,bVf] = svd(BOLD_Map(isub,:),'econ');
    
    [bUmixed,bDmixed,bVmixed] = svd(bDm'*bUm'*bUf*bDf,'econ');
    bSalience_M =bVm*bUmixed;
    bSalience_f = bVf*bVmixed;
    bDcorr = diag(bDmixed);
    
    
    % finding bs components that have the highest correlation with the 
    % original components and sort bs components based on their (pairwise) 
    % correlation value with original components
    
    % soritng bootsrap components to match the original components, 
    % (the bs comp with the highest correlation to original comp 1 is
    % assumed to be it's corresponding comp )
    %[indf,sgf] = sort_eigen_images(Salience_f,bSalience_f); %~ commented by Nasim 
    [indf,sgf] = sort_eigen_images(Salience_M,bSalience_M); % ~ added by Nasim 
    
    % adiing back the sign ( of the correlation values) to the components
    tm = bsxfun(@times,bSalience_M(:,indf),sgf);
    tf = bsxfun(@times,bSalience_f(:,indf),sgf);
    
    RSalience_M = RSalience_M + tm.^2;
    RSalience_f = RSalience_f + tf.^2;
    MSalience_M = MSalience_M + tm;
    MSalience_f = MSalience_f + tf;
    
    % for each component, find the region/channel that has the highest 
    %  weight (* error in the function code) 
    % but it is not used in the results that I want
    vtm = average_over_voxels(tm,num_regions);
    vtf = average_over_voxels(tf,num_regions);

    vRSalience_M = vRSalience_M + vtm.^2;
    vRSalience_f = vRSalience_f + vtf.^2;
    vMSalience_M = vMSalience_M + vtm;
    vMSalience_f = vMSalience_f + vtf;
    
    
end

% var(x) = E(x2) - E(x)^2
VSalience_M = RSalience_M/bs - (MSalience_M/bs).^2;
VSalience_f = RSalience_f/bs - (MSalience_f/bs).^2;
% bootstrap SE
VSalience_M = bs*VSalience_M/(bs-1);
VSalience_f = bs*VSalience_f/(bs-1);

ZSalience_M = Salience_M./sqrt(VSalience_M);
ZSalience_f = Salience_f./sqrt(VSalience_f);


vVSalience_M = vRSalience_M/bs - (vMSalience_M/bs).^2;
vVSalience_f = vRSalience_f/bs - (vMSalience_f/bs).^2;
vVSalience_M = bs*vVSalience_M/(bs-1);
vVSalience_f = bs*vVSalience_f/(bs-1);

vZSalience_M = average_over_voxels(Salience_M,num_regions)./sqrt(vVSalience_M);
vZSalience_f = average_over_voxels(Salience_f,num_regions)./sqrt(vVSalience_f);


function [ind,sg] = sort_eigen_images(Zm_full,Zm_in)
[numvoxels,numpcs_full] = size(Zm_full);
[numvoxels,numpcs_half] = size(Zm_in);

numpcs = min(numpcs_full,numpcs_half);

r_temp = corr(Zm_full,Zm_in);
rep_sign = sign(r_temp);
r_temp = abs(r_temp);

for i = 1:numpcs
    [iii,jjj] = find(r_temp == max(r_temp(:)));
    ind(iii) = jjj;
    sg(iii)  = sign(rep_sign(iii,jjj));
    
    r_temp(iii,:) = -1;
    r_temp(:,jjj) = -1;
end


function vZm = average_over_voxels(Zm,num_regions)

len = size(Zm,1)/num_regions;
if fix(len)~=len
    num_regions = size(Zm,1);
    len = 1;
end
vZm = zeros(num_regions,size(Zm,2));
for i = 1:num_regions
    tz = Zm((i-1)*len+1:i*len,:);
    [mx1,id1] = max(abs(tz));
    vZm(i,:) = tz(id1);
end





