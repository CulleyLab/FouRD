%% load input files [.tif] %%

warning off

% load folder where input files are stored
uiwait(msgbox('Load parent folder'));
superparent_d = uigetdir('');

recuseFlag = (questdlg('recurse subfolders?'));

matlab_folder = cd;
cd(superparent_d)
mkdir('output')
[parameters, listing_masks] = user_input(superparent_d);

% get list of subfolders
superFileList=dir;
dcount=0;
MeanAngleAll=[];
StdAngleAll=[];
MeanExccAll=[];
StdExccAll=[];
for fcount=3:numel(superFileList)
    if ((superFileList(fcount).isdir == 1) && (~strcmp(superFileList(fcount).name,'output')));
        fcount
        parent_d=strcat(superparent_d,'\',superFileList(fcount).name);
        cd(parent_d)
        listing = dir('*.tif');
        mkdir('output')
        cd(matlab_folder)
        [dirAnglemat,dirExccmat]=AFTbatchFunc(parent_d,listing,parameters);
        MeanAngleAll(:,:,fcount-2)=mean(dirAnglemat,3);
        StdAngleAll(:,:,fcount-2)=std(dirAnglemat,0,3);
        MeanExccAll(:,:,fcount-2)=mean(dirExccmat,3);
        StdExccAll(:,:,fcount-2)=std(dirExccmat,0,3);
    end
end



%% Section for functions
function [Anglemat,Exccentricitymat]=AFTbatchFunc(parent_d,listing,parameters) 
    n_files = length(listing);

    %limit range for testing
    %n_files=3;

    av_ordermat = zeros(n_files,1);

    for file_list = 1:n_files
    
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
        [av_ordermat(file_list,1),Anglemat(:,:,file_list),Exccentricitymat(:,:,file_list)] = AFT_function(file, directory, parameters);
    
    end

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
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'mean_angle_heatmap.tif'));
    close
    
    % std angle heat map
    figure;
    imagesc(rad2deg(std(Exccentricitymat,0,3)));
    %set(gca,'visible','off')
    caxis([0,90]);
    colormap(parula);
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'],'std_angle_heatmap.tif'));
    close
    
end