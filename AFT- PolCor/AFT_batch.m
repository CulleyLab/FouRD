%% load input files [.tif] %%

warning off

% load folder where input files are stored
uiwait(msgbox('Load parent folder'));
parent_d = uigetdir('');

matlab_folder = cd;
cd(parent_d)
listing = dir('*.tif');

% create output folder
mkdir('output')

cd(matlab_folder)

%% user input and set parameters %%

[parameters, listing_masks] = user_input(parent_d);

%% open one file at a time and perform analysis %%

n_files = length(listing);

%limit range for testing
%n_files=3;

av_ordermat = zeros(n_files,1);

% get info on file list to get image sizes
MaxFWidth = 0;
MaxFHeight = 0;
for file_list = 1:n_files
    frame_info = imfinfo(strcat(listing(file_list).folder,'\',listing(file_list).name));
    if (frame_info.Width > MaxFWidth)
        MaxFWidth = frame_info.Width;
    end
    if (frame_info.Height > MaxFHeight)
        MaxFHeight = frame_info.Height;
    end
end
%make matricies large enogh for bigest image
MaxFHeight=ceil(MaxFHeight/(parameters.winsize*parameters.overlap));
MaxFWidth=ceil(MaxFWidth/(parameters.winsize*parameters.overlap));
Anglemat = zeros(MaxFHeight,MaxFWidth,n_files);
Exccentricitymat = Anglemat;


tic
if (parameters.parproc==1)
    % predifine some variables for division to workers
    directory = listing(1).folder;

    parfor file_list = 1:n_files
        fprintf('Analysing file %d of %d',file_list,n_files)
        fprintf('\n')

        % file and directory name
        file = listing(file_list).name;

        % call function
        [atemp,btemp,ctemp]=AFT_function(file, directory, parameters);
        % make local matrix to store values in uper left corner
        % this is because of restricted inexing for sliced variables in a parfor loop
        bmaxi = size(btemp,1);
        bmaxj = size(btemp,2);
        Btemp = zeros(MaxFHeight,MaxFWidth);
        Ctemp = zeros(MaxFHeight,MaxFWidth);
        Btemp(1:bmaxi,1:bmaxj)=btemp;
        Ctemp(1:bmaxi,1:bmaxj)=ctemp;
        % update matricies with function data for current image
        av_ordermat(file_list,1)=atemp;
        Anglemat(:,:,file_list)=Btemp;
        Exccentricitymat(:,:,file_list)=Ctemp;
        
        %[av_ordermat(file_list,1),Anglemat(:,:,file_list),Exccentricitymat(:,:,file_list)] = AFT_function(file, directory, parameters);
    end
else
    for file_list = 1:n_files
        fprintf('Analysing file %d of %d',file_list,n_files)

        % file and directory name
        file = listing(file_list).name;
        directory = listing(file_list).folder;

        % file and directory name for mask (if local)
        if parameters.mask_method == 1
            file_mask = listing_masks(file_list).name;
            directory_mask = listing_masks(file_list).folder;
            parameters.mask_name = fullfile(directory_mask,file_mask);
        end
        
        % call function
        [atemp,btemp,ctemp]=AFT_function(file, directory, parameters);
        bmaxi = size(btemp,1);
        bmaxj = size(btemp,2);
        % update matricies with function data for current image
        av_ordermat(file_list,1)=atemp;
        Anglemat(1:bmaxi,1:bmaxj,file_list)=btemp;
        Exccentricitymat(1:bmaxi,1:bmaxj,file_list)=ctemp;
       
        %[av_ordermat(file_list,1),Anglemat(:,:,file_list),Exccentricitymat(:,:,file_list)] = AFT_function(file, directory, parameters);
     end
end
toc

% shrink Anglemat and excentricity mat to fit all non zero data
[nzi,nzj,nzk]=ind2sub([size(Anglemat,1) size(Anglemat,2) size(Anglemat,3)],find(Exccentricitymat>0));
nzi=max(nzi);
nzj=max(nzj);
nzk=max(nzk);
Anglemat=Anglemat(1:nzi,1:nzj,1:nzk);
Exccentricitymat=Exccentricitymat(1:nzi,1:nzj,1:nzk);


% save order parameter and maps
save(fullfile([parent_d '/output'], 'median_order_parameter.mat'), 'av_ordermat');
save(fullfile([parent_d '/output'], 'Anglemat.mat'), 'Anglemat');
save(fullfile([parent_d '/output'], 'Exccentricitymat.mat'), 'Exccentricitymat');

T = table(av_ordermat);
T.Properties.VariableNames = {'median_order_parameter'};
writetable(T,fullfile([parent_d '/output'], 'median_order_parameter.csv'))

% mean eccentricity heat map
    figure;
    imagesc(mean(Exccentricitymat,3));
    title('Mean Eccentricity');
    %set(gca,'visible','off');
    caxis([0,1]);
    colormap(parula);
    colorbar();
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'mean_excc_heatmap.tif'));
    close

%std eccentricity heat map
    figure;
    imagesc(std(Exccentricitymat,0,3));
    title('STD Eccentricity');
    %set(gca,'visible','off');
    caxis([0,1]);
    colormap(parula);
    colorbar();
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'std_excc_heatmap.tif'));
    close
    
% rel std eccentricity heat map
    figure;
    imagesc(std(Exccentricitymat,0,3)./mean(Exccentricitymat,3));
    title('STD/Mean Eccentricity');
    %set(gca,'visible','off');
    caxis([0,1]);
    colormap(parula);
    colorbar();
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'rel_mean_excc_heatmap.tif'));
    close
    
% mean angle heat map
    figure;
    imagesc(rad2deg(mean(Anglemat,3)));
    title('Mean Angle');
    %set(gca,'visible','off');
    caxis([0,180]);
    colormap(hsv);
    colorbar();
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'mean_angle_heatmap.tif'));
    close
    
% std angle heat map
    figure;
    imagesc(rad2deg(std(Anglemat,0,3)));
    title('STD Angle');
    %set(gca,'visible','off')
    caxis([0,90]);
    colormap(parula);
    colorbar();
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'std_angle_heatmap.tif'));
    close
    
%clear; clc