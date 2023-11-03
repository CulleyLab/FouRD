% This code is a modified version of the 'Alignment by Fourier Transform' package associated with the paper
% "A Workflow for Rapid Unbiased Quantification of Fibrillar Feature Alignment in Biological Images" Front. Comput. Sci 14 October 2021, Stefania Marcotti et.al.
% DOI: https://doi.org/10.3389/fcomp.2021.745831
% Original software avalible at https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform
% All rights and permissions belong to: Patrick Oakes poakes@gmail.com
% The modifications support intigration into the FouRD package submitied with the paper
% "Characterisation and correction of polarisation effects in fluorescently labelled fibres", Nandini Aggarwi et.al.
% submitted to the Journal of Microscopy. 3/11/2023

function [parameters, listing_masks] = user_input(parent_d)

% user input dialog
prompt = {'Window size [px]', ...
    'Window overlap [%]',...
    'Neighbourhood radius [vectors]'};
prompt_title = 'Parameters';
dims = [1 50];
definput = {'250','50','2'};
user_answer = inputdlg(prompt, prompt_title, dims, definput);

parameters_save.winsize_px = str2double(user_answer{1,1});
parameters_save.overlap_percentage = str2double(user_answer{2,1});
parameters_save.neighbourhood_radius = str2double(user_answer{3,1});

parameters.winsize = parameters_save.winsize_px;
parameters.overlap = 1 - parameters_save.overlap_percentage/100;
parameters.st = parameters_save.neighbourhood_radius;

% user answers on saving and filtering
answer_output = questdlg('Would you like to save output images?', ...
    'Output images', 'Yes', 'No', 'Yes');
switch answer_output
    case 'Yes'
        output_images = 1;
    case 'No'
        output_images = 0;
end

parameters_save.figures = output_images;
parameters.figures = parameters_save.figures;

% option for colorising vetors by eccentricity
answer_output = questdlg('Would you like to colorise alignment vectors by eccentricity?', ...
    'Colorise vectors', 'Yes', 'No', 'Yes');
switch answer_output
    case 'Yes'
        col_vectors = 1;
    case 'No'
        col_vectors = 0;
end

parameters_save.colvectors = col_vectors;
parameters.colvectors = parameters_save.colvectors;

% option for parralel processing
answer_parallel = questdlg('Use parallel processing? (Masks can not be used! and live results suppressed)', ...
    'Parallel processing','Yes', 'No', 'No');
switch answer_parallel
    case 'Yes'
        parallel_proc = 1;
    case 'No'
        parallel_proc = 0;
end

parameters_save.parproc = parallel_proc;
parameters.parproc = parameters_save.parproc;

% option for masks/filtering
if (parallel_proc == 0)
    answer_filter = questdlg('Would you like to apply local masking and/or filtering on the images?', ...
        'Filter images', 'Yes', 'No', 'Yes');
    switch answer_filter
        case 'Yes'
            filter_images = 1;
        case 'No'
            filter_images = 0;
    end
else
    filter_images = 0;
end

% initialise parameters if not going into if statement
parameters.checkpoint = 0;
parameters.mask_method = 0;
parameters.filter_blank = 0;
parameters.filter_isotropic = 0;
parameters.eccentricity_threshold = 0;
listing_masks = [];

if filter_images == 1
    prompt_filter = {'Masking method (local = 1, global = 0)', ...
        'Ignore blank spaces (yes = 1, no = 0)',...
        'Ignore isotropic regions (yes = 1, no = 0)'};
    prompt_filter_title = 'Filters';
    definput_filter = {'0','0','0'};
    user_answer_filter = inputdlg(prompt_filter, prompt_filter_title, dims, definput_filter);
    
    parameters_save.mask_method = str2double(user_answer_filter{1,1});
    parameters_save.filter_blank = str2double(user_answer_filter{2,1});
    parameters_save.filter_isotropic = str2double(user_answer_filter{3,1});
    
    parameters.mask_method = parameters_save.mask_method;
    parameters.filter_blank = parameters_save.filter_blank;
    parameters.filter_isotropic = parameters_save.filter_isotropic;
    
    % load masks if masking method == local
    if parameters.mask_method == 1
        
        uiwait(msgbox('Load folder containing masks'));
        masks_d = uigetdir('');
        
        matlab_folder = cd;
        cd(masks_d)
        listing_masks = dir('*.tif');
        
        cd(matlab_folder)
    end
    
    % if blank regions to be discarded
    if str2double(user_answer_filter{2,1}) == 1
        
        prompt_blank = {'Set mean intensity value of blank regions (0:black - 255:white)'};
        prompt_blank_title = 'Blank regions';
        definput_blank = {'0'};
        user_answer_blank = inputdlg(prompt_blank, prompt_blank_title, dims, definput_blank);
        
        parameters_save.filter_blank_value = str2double(user_answer_blank{1,1});
        parameters.checkpoint = parameters_save.filter_blank_value / 255; % threshold mean in each window to do the calculation
        
    end
    
    % if isotropic regions to be discarded
    if str2double(user_answer_filter{3,1}) == 1
        
        prompt_isotropic = {'Set threshold for FFT eccentricity (0:isotropic - 1:highly oriented)'};
        prompt_isotropic_title = 'Isotropic regions';
        definput_isotropic = {'0'};
        user_answer_isotropic = inputdlg(prompt_isotropic, prompt_isotropic_title, dims, definput_isotropic);
        
        parameters_save.filter_isotropic_value = str2double(user_answer_isotropic{1,1});
        parameters.eccentricity_threshold = parameters_save.filter_isotropic_value; % threshold FFT eccentricity in each window to do the calculation
        
    end
    
end

% save parameters in [output] folder
save(fullfile([parent_d '/output'], 'parameters.mat'), 'parameters_save');

T = struct2table(parameters_save);
writetable(T,fullfile([parent_d '/output'], 'parameters.txt'))

end