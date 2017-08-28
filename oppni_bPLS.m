
%% check if all input variables have been provided %%

error_quit = 0 ;

if ( exist('OPPNI_dir') == false ) || ( exist(OPPNI_dir,'dir') == false )
	error('Must enter a valid OPPNI directory.') ;
	error_quit = error_quit + 1 ;
end

if ( exist('BEHAV_data') == false ) || ( exist(BEHAV_data,'dir') == false )
	error('Must enter a valid behaviour file.') ;
	error_quit = error_quit + 1 ;
end

if ( exist('OUTPUT_dir') == false ) 
	error('Must provide an output directory.')
	error_quit = error_quit + 1 ;
end

if (exist('prefix') == false) && (error_quit == 0)
	prefix = 'output' ;
	error('No prefix provided. Defaulting to "output". ')
end

if (exist('pipeline') == false) || (pipeline < 1) || (pipeline > 3)
	error( sprintf(['You must provide a valid pipeline value \n '...
									'		1 = CON \n' ...
									'		2 = FIX \n' ...
									'		3 = IND \n']) )
	error_quit = error_quit + 1 ;
end

if (exist('var_norm') == false) || (var_norm < 0) || (var_norm > 2)
	error( sprintf(['You must provide a valid pipeline value. \n '...
									'		0 = no normalization, directly use input values \n' ...
									'		1 = (column wise) mean centring  of X and Y \n' ...
									'		2 = zscore X and Y  \n']) )
	error_quit = error_quit + 1 ;
end

% check if any errors have occured
if error_quit > 0
	exit 
end

%% set output filename %%

prefix = [prefix, '_bPLS.mat'] ;

%% import behavioural data %%

mats.Y = dlmread(BEHAV_data) ;
if prod( size(mats.Y) ) == 0
	error('Behaviour file contains no data.')
	exit
end

%% index all cSPM files %%

SPM_dir = fullfile(OPPNI_dir, 'optimization_results', 'spms') ;
if ( exist(OPPNI_dir,'dir') == false )
	error('The spm subdirectory does not exist in the OPPNI output directory provided.')
	exit
end

cSPMs = dir( fullfile(SPM_dir, '*_sNorm.nii') ) ;
cSPMs = {cSPMs(:).name}' ;

if length(cSPMs) == 0
	error('No spm files exist.')
	exit
end

if size(mats.Y,1) ~= size(cSPMs,1)
	error('Behaviour file must have the same number of rows as there are cSPM files.')
	exit
end

%% loop to import functional data %%

count = 0 ;
for cSPM = cSPMs

	count = count + 1 ;

	cSPM = fullfile(SPM_dir, cSPM) ;
	cSPM = load_nii(cSPM) ;
	cSPM = cSPM.nii ;
	cSPM = cSPM(:,:,:, pipeline) ;

	dims = prod( size(cSPM) ) ;
	dims = [ 1, dims ] ;

	cSPM = reshape( cSPM, dims ); 

	v_cSPM = reshape( cSPM, dims ) ;
	
	mats.X(count, :) = v_cSPM ;

end

%% running behavioural PLS analysis %%

[plsResults.avg_ZSalience_X, plsResults.avg_ZSalience_Y, plsResults.pred_scores_X , plsResults.pred_scores_Y, plsResults.pls_out] = pls_nasim(mats.X, mats.Y, var_norm) ;

%% saving results %%

OUTPUT_path = fullfile(OUTPUT_dir, prefix) ;

save(OUTPUT_path, plsResults) ;

exit
