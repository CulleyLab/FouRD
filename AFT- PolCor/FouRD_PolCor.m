%% Polarisation correction utilities
%Corrects for polarisation bias of aligned samples in microscopy images when dye labels are
%not fully free to rotate
%Requires a calibration sample (2d image or Z stack) that has no intrinsic
%net alignment. The anisotropy of each ring in Fourier space is measured to
%get calibrtion data. The inverse of this spectrum can then be appled to an
%aligned sample and the effect of the polarisation bias removed to leave
%just the true alignment. 

%% PRogram to creat user interface to lanch programs related to polarisation correction
%GUI parameters
figSize = [200 280];
figMargin = 20;
panelMargin = 10;
buttonWidth = 180;
buttonHeight = 25;

%  Initialize and construct the GUI.
guifig = figure('Visible','on','Position',[100,100,figSize],...
    'MenuBar','None','ToolBar','none',...
    'DefaultUIPanelUnits','pixels','DefaultAxesUnits','pixels',...
    'DefaultUIPanelFontName','SegoeUI','DefaultUIPanelFontSize',10,...
    'DefaultUIControlFontName','SegoeUI','DefaultUIControlFontSize',9);
guifig.Name='PolCorUtils';

% add text
textlabel = uicontrol(guifig,'Style','text','String','Utilities to measure and correct for polarisation bias in aligned samples.');
textlabel.Position = [panelMargin, 220,buttonWidth,50];

% add buttons
buttonCalib = uicontrol(guifig,'Style','pushbutton','String','Make Calibration','position',[panelMargin,100,buttonWidth,buttonHeight]);
buttonCorrect = uicontrol(guifig,'Style','pushbutton','String','Apply correction','position',[panelMargin,70,buttonWidth,buttonHeight]);
buttonImageSpec = uicontrol(guifig,'Style','pushbutton','String','Spectrum from Image','position',[panelMargin,40,buttonWidth,buttonHeight]);
buttonCalibSpec = uicontrol(guifig,'Style','pushbutton','String','Spectrum from Calibration','position',[panelMargin,10,buttonWidth,buttonHeight]);
buttonBatchAFT = uicontrol(guifig,'Style','pushbutton','String','Batch AFT','position',[panelMargin,160,buttonWidth,buttonHeight]);
buttonSuperAFT = uicontrol(guifig,'Style','pushbutton','String','Recursive batch AFT','position',[panelMargin,130,buttonWidth,buttonHeight]);
buttonAniSpec = uicontrol(guifig,'Style','pushbutton','String','Anisotropy from Image pair','position',[panelMargin,190,buttonWidth,buttonHeight]);

% add callbacks - starts named function when button pressed
buttonCalib.Callback = @startCalib;
buttonCorrect.Callback = @startCorrect;
buttonImageSpec.Callback = @startImageSpec;
buttonCalibSpec.Callback = @startCalibSpec;
buttonBatchAFT.Callback = @startBatchAFT;
buttonSuperAFT.Callback = @startSuperAFT;
buttonAniSpec.Callback = @startAniSpec;

% callback functions
% Each function gets parent handle, disables bottons to prevent
% simultanious function call. displays operation in text label calls the
% program associated with the operation. When finish reenables buttons

% callback function for calibration.
function startCalib(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Calculating a caibration matrix from stack.')
    finished=doGetCalibMat();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end

% callback for correction
function startCorrect(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Correcting a stack from a previous calibration.')
    finished=doCorFromCalib();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end

%callback for anisotropy spectrum from image.
function startImageSpec(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Calculating Fourier anisotropy spectrom of an image.')
    finished=doSpecFromFrame();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end

% callback for anisotropy spectrum from previous calibration data
function startCalibSpec(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Calculating Fourier anisotropy spectrom of a calibration matrix.')
    finished=doSpecFromCalib();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end

% callback to start AFT analysis of a directory
function startBatchAFT(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Perfoming AFT analysisof a stack')
    AFT_batch();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end

%callback to start AFT recursing sub directories and gathering statistics
function startSuperAFT(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Performing AFT analysis recusing through sub directories.')
    AFT_Superbatch();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end

% callback to get anisotropy spectrum from polarisation resolved image pair
function startAniSpec(src,event)
    fprintf('\n')
    Pfig=src.Parent;
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','off');
    end
    
    set(Pfig.Children(8),'String','Calculating Anisotropy from image pair.')
    finished=doAniFromPair();
    set(Pfig.Children(8),'String','Utilities to measure and correct for polarisation bias in aligned samples.')
    
    for encount = 1:7
        set(Pfig.Children(encount),'Enable','on');
    end
end


%% These functions perform the various operations for each callback. each operates independantly, any common code is in the 'required functions' section at the end
%each are called from the callbacks above and return a flag if it compleated succesfully or aborted.

%% program to get calibration anisotropy from a z stack
function finished = doGetCalibMat()
% load folder where calibration files are stored
uiwait(msgbox('Load calibration folder'));
parent_d = uigetdir('');
if (parent_d==0)    % abort if cancelled
    disp('Aborting calibration')
    finished=false;
    return
end
matlab_folder = cd;
cd(parent_d)
listing = dir('*.tif');

% user input dialog for calibration parameters
prompt = {'Max Pixels X', 'Max Pixels Y', ...                   % crops images to this size (usefull if images are not all the same size from merging sub scans)
    'Number of tiles X', ...                                    % Number of subdivitions of calibration image in X dimention if tiling
    'Number of tiles Y',...                                     % Number of subdivisions of calibration image in Y direction if tiling   
    'Start ring'};                                              % Spatial frequency to start the corection from (below this is uncorrected)
prompt_title = 'Calibration Parameters';
dims = [1 50];
definput = {'1024','1024','7','7','2'};
user_answer = inputdlg(prompt, prompt_title, dims, definput);
if (isempty(user_answer) == 1)    % abort if cancelled
    disp('Aborting calibration')
    finished=false;
    return
end

% get and save parameter values
parameters.numPixelsX = str2double(user_answer{1,1});
parameters.numPixelsY = str2double(user_answer{2,1});
parameters.numTilesX = str2double(user_answer{3,1});
parameters.numTilesY = str2double(user_answer{4,1});
parameters.startRing = str2double(user_answer{5,1});

numpixelsX=parameters.numPixelsX;
numpixelsY=parameters.numPixelsY;
numtilesX=parameters.numTilesX;
numtilesY=parameters.numTilesY;

% initiate matricies to store moments
numfiles = length(listing);
calibP00=zeros(numtilesY,numtilesX,10000);
calibP20=zeros(numtilesY,numtilesX,10000);
calibP22=zeros(numtilesY,numtilesX,10000);

% loop over images in directory list (usually all the images of a z stack but can be a single image)
for filecount = 1:numfiles
    % frame progress indicator
    disp(strcat('Calibrating frame: ',num2str(filecount),' out of: ',num2str(numfiles)))
    
    % file and directory name
    file = listing(filecount).name;
    directory = listing(filecount).folder;
    
    % load image
    calibframe = imread(fullfile(directory, file));
    calibframe = im2double(calibframe);
    
    %clip image if exceeds max size (for irregular sized images)
    if (size(calibframe,1) > numpixelsY)
        calibframe=calibframe(1:numpixelsY,:);
    end
    if (size(calibframe,2) > numpixelsX)
        calibframe=calibframe(:,1:numpixelsX);
    end
    
    %format and offset of tile progress indicator
    fprintf('                                    \n')
    formspec='Looping over tiles: Y = %3d, X = %3d';
    backmsg=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
    %loop over tiles geting start and end pixel indicies
    for ycount = 1:numtilesY
        for xcount = 1:numtilesX
            %tile progress indicator
            message=sprintf(formspec,ycount,xcount);
            disp([backmsg message])
            
            % calculate start and end indicies for tile
            startX=((xcount-1)*floor(size(calibframe,2)/numtilesX))+1;
            startY=((ycount-1)*floor(size(calibframe,1)/numtilesY))+1;
            
            if (xcount==numtilesX && ycount==numtilesY) % if last x&y tile
                endX=size(calibframe,2);
                endY=size(calibframe,1);
            elseif (xcount==numtilesX) % if last x tile
                endX=size(calibframe,2);
                endY=((ycount-0)*floor(size(calibframe,1)/numtilesY))+1;
            elseif (ycount==numtilesY) % if last y tile
                endY=size(calibframe,1);
                endX=((xcount-0)*floor(size(calibframe,2)/numtilesX))+1;
            else % if not edge tile
                endX=((xcount-0)*floor(size(calibframe,2)/numtilesX))+1;
                endY=((ycount-0)*floor(size(calibframe,1)/numtilesY))+1;
            end
            
            % get anisotropy of tile and update totals (calls common spectrum function in required functions section bellow)
            calibpatch=calibframe(startY:endY,startX:endX);
            [rawPF,rawAF,rawCor,rawP00,rawP20,rawP22]=getImgSpectrum(calibpatch,false);
            calibP00(ycount,xcount,1:numel(rawP00))=calibP00(ycount,xcount,1:numel(rawP00))+reshape(rawP00,[1 1 numel(rawP00)]);
            calibP20(ycount,xcount,1:numel(rawP20))=calibP20(ycount,xcount,1:numel(rawP20))+reshape(rawP20,[1 1 numel(rawP20)]);
            calibP22(ycount,xcount,1:numel(rawP22))=calibP22(ycount,xcount,1:numel(rawP22))+reshape(rawP22,[1 1 numel(rawP22)]);
        end
    end
    
    % move cursor back to overwrite previous progress indicator
    disp(backmsg)
    fprintf('\b')
end

%create calibration matrix and save
clear calibmat
calibmat(:,:,:,1)=calibP00/numfiles;
calibmat(:,:,:,2)=calibP20/numfiles;
calibmat(:,:,:,3)=calibP22/numfiles;
calibmat=calibmat(:,:,1:numel(rawP00),:);
matrixname=strcat('calibrationmat_X',num2str(size(calibmat,2)),'Y',num2str(size(calibmat,1)),'.mat');
save(matrixname,'calibmat');
cd (matlab_folder);

%report successfull compleation
fprintf('\nDone geting calibration data from stack.\n')
finished = true;
end

%% Programe to correct a stack from a calibration file.
function finished = doCorFromCalib()
% get calibration matrix and parameters from previous calibration
uiwait(msgbox('Load calibration data'));
[calibfile,calibfolder] = uigetfile('*.mat');
if (calibfile==0)    % abort if cancelled
    disp('Aborting correction')
    finished=false;
    return
end
matlab_folder = cd;
cd(calibfolder);
clear calibmat
load(calibfile);
cd(matlab_folder);
parameters.calibfolder=calibfolder;

% get stack for corection.
uiwait(msgbox('select stack folder for correction'));
parent_d = uigetdir('');
if (parent_d==0)    % abort if cancelled
    disp('Aborting correction')
    finished=false;
    return
end
cd(parent_d)
listing = dir('*.tif');
% get parameters from previous calibration
parameters.numTilesX = size(calibmat,2);
parameters.numTilesY = size(calibmat,1);
numtilesX=parameters.numTilesX;
numtilesY=parameters.numTilesY;

% user input dialog for corection params
prompt = {'Max Pixels X', 'Max Pixels Y'};              % crops images to this size (Shold be the same as calibration data)
prompt_title = 'Analysis Parameters';
dims = [1 50];
definput = {'1024','1024'};
user_answer = inputdlg(prompt, prompt_title, dims, definput);
if (isempty(user_answer)==1)    % abort if cancelled
    disp('Aborting correction')
    finished=false;
    return
end
% get size parameters from user dialog
parameters.numPixelsX = str2double(user_answer{1,1});
parameters.numPixelsY = str2double(user_answer{2,1});
numpixelsX=parameters.numPixelsX;
numpixelsY=parameters.numPixelsY;

% user input dialog
prompt = {strcat('Tiles X=: ',num2str(numtilesX),' ,Tiles Y=: ',num2str(numtilesY),' ,Start ring :'),'Rings to smooth :'}; % start frequency for correction & degree of smoothing of calibration spectrum.
prompt_title = 'Corection Parameters';
dims = [1 50];
definput = {'2','1'};
user_answer = inputdlg(prompt, prompt_title, dims, definput);
if (isempty(user_answer)==1)    % abort if cancelled
    disp('Aborting correction')
    finished=false;
    return
end
% get convertion parameters from user dialog
parameters.startRing = str2double(user_answer{1,1});
parameters.smoothFac = str2double(user_answer{2,1});
startring=parameters.startRing;
smoothfac=parameters.smoothFac;

%make target dir
correct_d=strcat(parent_d,'_CorrX',num2str(numtilesX),'Y',num2str(numtilesY),'SR',num2str(startring),'SF',num2str(smoothfac));
mkdir(correct_d);
cd(correct_d);

% initials matricies for moments
numfiles = length(listing);
calibP00=calibmat(:,:,:,1);
calibP20=calibmat(:,:,:,2);
calibP22=calibmat(:,:,:,3);

% loop over all files in the directory list (Usually a Z stack but can contain all the measuments for a common calibration)
for filecount = 1:numfiles
    %frame progress indicator
    disp(strcat('Correcting frame: ',num2str(filecount),' out of: ',num2str(numfiles)))
    
    % file and directory name
    file = listing(filecount).name;
    directory = listing(filecount).folder;
    % load image
    imageframe = imread(fullfile(directory, file));
    
    %clip image if exceeds max size (for irregular sized images)
    if (size(imageframe,1) > numpixelsY)
        imageframe=imageframe(1:numpixelsY,:);
    end
    if (size(imageframe,2) > numpixelsX)
        imageframe=imageframe(:,1:numpixelsX);
    end
    
    %format and offset of tile progress indicator
    fprintf('                                    \n')
    formspec='Looping over tiles: Y = %3d, X = %3d';
    backmsg=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');

    %loop over tiles getinfg start and end pixel indicies
    corframe=imageframe*0;
    for ycount = 1:numtilesY
        for xcount = 1:numtilesX
            %tile progress indicator
            message=sprintf(formspec,ycount,xcount);
            disp([backmsg message])
            
            % calculate start and end indicies of current tile as above
            startX=((xcount-1)*floor(size(imageframe,2)/numtilesX))+1;
            startY=((ycount-1)*floor(size(imageframe,1)/numtilesY))+1;
            if (xcount==numtilesX && ycount==numtilesY)
                endX=size(imageframe,2);
                endY=size(imageframe,1);
            elseif (xcount==numtilesX)
                endX=size(imageframe,2);
                endY=((ycount-0)*floor(size(imageframe,1)/numtilesY))+1;
            elseif (ycount==numtilesY)
                endY=size(imageframe,1);
                endX=((xcount-0)*floor(size(imageframe,2)/numtilesX))+1;
            else
                endX=((xcount-0)*floor(size(imageframe,2)/numtilesX))+1;
                endY=((ycount-0)*floor(size(imageframe,1)/numtilesY))+1;
            end
            
            %get corection data and corect patch
            rawP00=calibP00(ycount,xcount,:);
            rawP20=calibP20(ycount,xcount,:);
            rawP22=calibP22(ycount,xcount,:);
            imagepatch=imageframe(startY:endY,startX:endX);
            corpatch=corPolImgFromSpec(imagepatch,rot90(smooth(rawP00,smoothfac),1),rot90(smooth(rawP20,smoothfac),1),rot90(smooth(rawP22,smoothfac),1),startring,false);
            corframe(startY:endY,startX:endX)=corpatch;
        end
    end
    
    % move cursor back to overwrite previous progress indicator
    disp(backmsg)
    fprintf('\b')
    
    % save corected image
    simpic(:,:,1)=uint16(corframe);
    imwrite(simpic(:,:,1),strcat(correct_d,'\',file),'WriteMode','overwrite','Compression', 'none');
end

%save parameters used
save('parameters.mat','parameters')
cd(matlab_folder)

%report successfull compleation
fprintf('\nDone corecting stack from calibration data.\n')
finished=true;
end

%% program to get Forier ring anisotropy from image frame
function finished = doSpecFromFrame()
% load image
uiwait(msgbox('Load image file'));
[imagefile,imagefolder] = uigetfile('*.tif');
if (imagefile==0)   % abort if cancelled
    disp('Aborting calculation')
    finished=false;
    return
end
imageframe = imread(fullfile(imagefolder, imagefile));

%calculate spectrum (calls common spectrum function in required functions section bellow)
disp('Calculating anisotropy spectrum.')
[rawPF,rawAF,rawCor,rawP00,rawP20,rawP22]=getImgSpectrum(imageframe,true);

%plot graph
figa=figure;
plot(rawP20./rawP00,'c','LineWidth',1)
hold on
plot(rawP22./rawP00,'m','LineWidth',1')
plot(smooth(rawP20./rawP00,10),'b','LineWidth',2)
plot(smooth(rawP22./rawP00,10),'r','LineWidth',2)
legend('Y20','Y22','Y20smooth','Y22smooth')
xlabel('spacial frequency')
ylabel('value')
title('Anisotropy spectrum')

%report successful compleation
fprintf('\nDone\n')
finished=true;
end

%% program to get Furrier ring anisotropy from a polarisation resolved image pair 
function finished = doAniFromPair()
% load parallel polarised image
uiwait(msgbox('Load parallel polarised image file'));
[imagefile,imagefolder] = uigetfile('*.tif');
if (imagefile==0)   %abort if cancelled
    disp('Aborting calculation')
    finished=false;
    return
end
parframe = imread(fullfile(imagefolder, imagefile));

% load perpendicular polarised image
uiwait(msgbox('Load perpendicular polarised image file'));
[imagefile,imagefolder] = uigetfile('*.tif');
if (imagefile==0)   %abort if cancelled
    disp('Aborting calculation')
    finished=false;
    return
end
perpframe = imread(fullfile(imagefolder, imagefile)); 

%do gfactor correction
prompt = {'G factor to multiply perpendicular image by?'};      % factor to multiply second channel by to get equal detection efficiency
prompt_title = 'Gfactor';
dims = [1 50];
definput = {'1'};
user_answer = inputdlg(prompt, prompt_title, dims, definput);
if (isempty(user_answer)==1)    % about if cancelled
    disp('Aborting correction')
    finished=false;
    return
end
perpframe=perpframe*str2double(user_answer);

%calculate spectrum ( uses common PairSpectrum function from required functions section bellow)
disp('Calculating pair anisotropy spectrum.')
[rawPF,rawAF,rawCor,rawP00,rawP20,rawP22]=getPairSpectrum(parframe,perpframe,true);

%plot graph
figa=figure;
plot(rawP20./rawP00,'c','LineWidth',1)
hold on
plot(rawP22./rawP00,'m','LineWidth',1')
plot(smooth(rawP20./rawP00,10),'b','LineWidth',2)
plot(smooth(rawP22./rawP00,10),'r','LineWidth',2)
legend('Y20','Y22','Y20smooth','Y22smooth')
xlabel('spacial frequency')
ylabel('value')
title('Anisotropy spectrum')

%report sucessful compleation
fprintf('\nDone\n')
finished=true;
end

%% program to get Furrier ring anisotropy from calibration matrix
function finished = doSpecFromCalib()
% get calibration matrix
uiwait(msgbox('Load calibration data'));
[calibfile,calibfolder] = uigetfile('*.mat');
if (calibfile==0)   %abort if cancelled
    disp('Aborting calculation')
    finished=false;
    return
end
matlab_folder = cd;
cd(calibfolder);
clear calibmat
load(calibfile);
cd(matlab_folder);
% get parameters from calibration
parameters.calibfolder=calibfolder;
parameters.numTilesX = size(calibmat,2);
parameters.numTilesY = size(calibmat,1);
numtilesX=parameters.numTilesX;
numtilesY=parameters.numTilesY;

% user input dialog
prompt = {strcat('Tiles X=: ',num2str(numtilesX),' ,Tiles Y=: ',num2str(numtilesY),' ,Rings to smooth :'),...   % amount of smoothing to spectrum & display a seperat spectrum for each tile or just qhole image
    'Tile ? :'};
prompt_title = 'Display Parameters';
dims = [1 50];
definput = {'3','N'};
user_answer = inputdlg(prompt, prompt_title, dims, definput);
if (isempty(user_answer)==1)    %abort if cancelled
    disp('Aborting calculation')
    finished=false;
    return
end
% get parameters from dialoge
parameters.smoothFac = str2double(user_answer{1,1});
parameters.tileSpectrum=user_answer{2,1};
smoothfac=parameters.smoothFac;
if (parameters.tileSpectrum=='Y')
    tilespectrum=true;
else
    tilespectrum=false;
end

%loop over tiles geting spectrum
disp('Calculating anisotropy spectrum.')
if (tilespectrum==true)     % if displaying seperate graphs    
    %loop over each tile and produce seperate spectrum for each tile
    for ycount = 1:numtilesY
        for xcount = 1:numtilesX
            %get anisotropy spectrum from calibration matrix
            rawP00=squeeze(calibmat(ycount,xcount,:,1));
            rawP20=squeeze(calibmat(ycount,xcount,:,2));
            rawP22=squeeze(calibmat(ycount,xcount,:,3));
            
            %plot graph
            figure
            plot(rawP20./rawP00,'c','LineWidth',1)
            hold on
            plot(rawP22./rawP00,'m','LineWidth',1')
            plot(smooth(rawP20./rawP00,smoothfac),'b','LineWidth',2)
            plot(smooth(rawP22./rawP00,smoothfac),'r','LineWidth',2)
            legend('Y20','Y22','Y20smooth','Y22smooth')
            xlabel('spacial frequency')
            ylabel('value')
            title(strcat('Anisotropy spectrum X=: ',num2str(xcount),' ,Y=: ',num2str(ycount)));
        end
    end
else    % just a single spectrum for the combined image
    % get anisotropy spectrum averaged over each tile in calibration
    rawP00=squeeze(mean(mean(calibmat(:,:,:,1),1),2));
    rawP20=squeeze(mean(mean(calibmat(:,:,:,2),1),2));
    rawP22=squeeze(mean(mean(calibmat(:,:,:,3),1),2));
    
    % plot the graph
    figure
    plot(rawP20./rawP00,'c','LineWidth',1)
    hold on
    plot(rawP22./rawP00,'m','LineWidth',1')
    plot(smooth(rawP20./rawP00,smoothfac),'b','LineWidth',2)
    plot(smooth(rawP22./rawP00,smoothfac),'r','LineWidth',2)
    legend('Y20','Y22','Y20smooth','Y22smooth')
    xlabel('spacial frequency')
    ylabel('value')
    title(strcat('Anisotropy spectrum tile average'));
end

%report succesful compleation
fprintf('\nDone producing tiled anisotropy spectrum.\n')
finished = true;
end

%% programe to do tiled correction of polarisation bias in an image frame
function finished = doCorTileFrame()
numtilesX=7;
numtilesY=7;
smoothfac=1;
startring=5;
cordata1=data1*0;

%loop over each tile
for ycount = 1:numtilesY
    for xcount = 1:numtilesX
        % get start and end indicies for each tile
        startX=((xcount-1)*floor(size(data1,2)/numtilesX))+1;
        startY=((ycount-1)*floor(size(data1,1)/numtilesY))+1;
        if (xcount==numtilesX && ycount==numtilesY)
            endX=size(data1,2);
            endY=size(data1,1);
        elseif (xcount==numtilesX)
            endX=size(data1,2);
            endY=((ycount-0)*floor(size(data1,1)/numtilesY))+1;
        elseif (ycount==numtilesY)
            endY=size(data1,1);
            endX=((xcount-0)*floor(size(data1,2)/numtilesX))+1;
        else
            endX=((xcount-0)*floor(size(data1,2)/numtilesX))+1;
            endY=((ycount-0)*floor(size(data1,1)/numtilesY))+1;
        end
        
        % get image patch of tile position 
        patch1=data1(startY:endY,startX:endX);
        
        % get spectrum from csalibration and aplly corection to patch (uses spectrum and corection fuctions from required functions section bellow)
        [rawPF,rawAF,rawCor,rawP00,rawP20,rawP22]=getImgSpectrum(patch1,true);
        corpatch1=corPolImgFromSpec(patch1,rot90(smooth(rawP00,smoothfac),1),rot90(smooth(rawP20,smoothfac),1),rot90(smooth(rawP22,smoothfac),1),startring,true);
        [corPF,corAF,corCor,corP00,corP20,corP22]=getImgSpectrum(corpatch1,true);
        
        %figure
        %plot(smooth(rawP20./rawP00,10),'b','LineWidth',1)
        %hold on
        %plot(smooth(rawP22./rawP00,10),'c','LineWidth',1)
        %plot(smooth(corP20./corP00,10),'r','LineWidth',1)
        %plot(smooth(corP22./corP00,10),'m','LineWidth',1)
        
        %add result to corected image
        cordata1(startY:endY,startX:endX)=corpatch1;
    end
end

%report successful compleation
fprintf('\nDone\n')
finished=true;
end

%% Section for required functions
% this section contains the common code to the above programs. These
% calculate the actual spectrum from an image patch or apply the corection
% to a patch given a set of anisotropy moments. Input args. are the images,
% the moments, any parameters and a flag indicating wether to update
% progress to the screen

% Function to measure the furier ring anisotropy of an image and return
% anissotropy spectum
function [polfac,angfac,corImg,P00,P20,P22]=getImgSpectrum(rawImg,progFlag)
    if (progFlag==true)
        fprintf('                         \n')
        formspec='Looping over rings: %3d %%';
        backmsg=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    end
    
    % Fourier space image
    FTmat=abs(fftshift(fft2(rawImg)));
    
    %generate angle and radius matricies
    [nx,ny]=meshgrid([-floor(size(FTmat,2)/2):1:ceil(size(FTmat,2)/2)-1],[-floor(size(FTmat,1)/2):1:ceil(size(FTmat,1)/2)-1]);
    if (size(FTmat,1)>size(FTmat,2))
        ny=ny*size(FTmat,1)/size(FTmat,1);
    end
    if (size(FTmat,2)>size(FTmat,1))
        nx=nx*size(FTmat,1)/size(FTmat,2);
    end
    Rmat=((nx.^2)+(ny.^2)).^0.5;
    Angmat=atan2d(ny,nx)+0.01;
    P00=[];
    P20=[];
    P22=[];
    N00=[];
    
    % loop over all radii (spatial frequencies)
    maxcount=min(floor(size(FTmat,2)/2),floor(size(FTmat,1)/2));
    maxcount=floor(maxcount*1.4*1);
    for icount= 1:maxcount-0
        temp=[];
     
        [i,j,v]=find(((Rmat>=icount) & (Rmat<(icount+1))).*FTmat);
        [i,j,ang]=find(((Rmat>=icount) & (Rmat<(icount+1))).*Angmat);
        N00(icount)=numel(i);
        ang=ang(1:numel(v));
        
        %calculate moments
        P00(icount)=sum(v);
        P00(icount)=P00(icount)/(2*3.142);
        P20(icount)=sum(v.*cosd(2*ang));
        P20(icount)=P20(icount)/3.142;
        P22(icount)=sum(v.*sind(2*ang));
        P22(icount)=P22(icount)/3.142;
        
        %progress indication
        if ((rem(icount*100/maxcount,10) < (100/maxcount)) && (progFlag==true))
            p=round(icount*100/maxcount);
            message=sprintf(formspec,p);
            disp([backmsg message])
        end
    end
    corImg=0;
    
    % calculate polarisation factor and angle
    temp=sqrt(((P20./P00).^2)+((P22./P00).^2));
    polfac=mean(temp(round(numel(temp)/8):round(numel(temp)/4)));
    temp=0.5*atan2d(P22,P20);
    angfac=mean(temp(round(numel(temp)/8):round(numel(temp)/4)));
    
    %progress
    if (progFlag==true)
        disp(backmsg)
        fprintf('\b')
    end
end


%% Function to measure the furier ring anisotropy from a pair of images and return
%anissotropy spectum
function [polfac,angfac,corImg,P00,P20,P22]=getPairSpectrum(rawImgA,rawImgB,progFlag)
    if (progFlag==true)
        fprintf('                         \n')
        formspec='Looping over rings: %3d %%';
        backmsg=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    end
    
    % Fourier space anisotropy image
    %FTmat=abs(fftshift(fft2(rawImg)));
    FTmat=abs(fft2(rawImgA))-abs(fft2(rawImgB));
    FTmat=FTmat./(abs(fft2(rawImgA))+abs(fft2(rawImgB)));
    FTmat=fftshift(FTmat+1);
    
    %generate angle and radii matricies
    [nx,ny]=meshgrid([-floor(size(FTmat,2)/2):1:ceil(size(FTmat,2)/2)-1],[-floor(size(FTmat,1)/2):1:ceil(size(FTmat,1)/2)-1]);
    if (size(FTmat,1)>size(FTmat,2))
        ny=ny*size(FTmat,1)/size(FTmat,1);
    end
    if (size(FTmat,2)>size(FTmat,1))
        nx=nx*size(FTmat,1)/size(FTmat,2);
    end
    Rmat=((nx.^2)+(ny.^2)).^0.5;
    Angmat=atan2d(ny,nx)+0.01;
    P00=[];
    P20=[];
    P22=[];
    N00=[];
    
    %Loop over al radii (spatial frequencies)
    maxcount=min(floor(size(FTmat,2)/2),floor(size(FTmat,1)/2));
    maxcount=floor(maxcount*1.4*1);
    for icount= 1:maxcount-0
        temp=[];
     
        [i,j,v]=find(((Rmat>=icount) & (Rmat<(icount+1))).*FTmat);
        [i,j,ang]=find(((Rmat>=icount) & (Rmat<(icount+1))).*Angmat);
        N00(icount)=numel(i);
        ang=ang(1:numel(v));
        
        % calculate moments
        P00(icount)=sum((v));
        P00(icount)=P00(icount)/(2*3.142);
        
        P20(icount)=sum(v.*cosd(2*ang));
        P20(icount)=P20(icount)/3.142;
        P22(icount)=sum(v.*sind(2*ang));
        P22(icount)=P22(icount)/3.142;
        
        %progress indication
        if ((rem(icount*100/maxcount,10) < (100/maxcount)) && (progFlag==true))
            p=round(icount*100/maxcount);
            message=sprintf(formspec,p);
            disp([backmsg message])
        end
    end
    corImg=0;
    
    % calculate polarisation factor and angle
    temp=sqrt(((P20./P00).^2)+((P22./P00).^2));
    polfac=mean(temp(round(numel(temp)/8):round(numel(temp)/4)));
    temp=0.5*atan2d(P22,P20);
    angfac=mean(temp(round(numel(temp)/8):round(numel(temp)/4)));
    
    %progress
    if (progFlag==true)
        disp(backmsg)
        fprintf('\b')
    end
end

%% Function to apply corection to image given Anisotropy spectrum
function [corImg] = corPolImgFromSpec(polImg,CP00,CP20,CP22,startFR,progFlag)
    if (progFlag==true)
        fprintf('                         \n')
        formspec='Looping over rings: %3d %%';
        backmsg=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    end
    
    %Fourier images
    %FTmat=abs(fftshift(fft2(polImg)));
    imgFT=fftshift(fft2(polImg));
    FTmat=abs(imgFT);
    %generate angle and radii matricies and blank inverse operator
    [nx,ny]=meshgrid([-floor(size(FTmat,2)/2):1:ceil(size(FTmat,2)/2)-1],[-floor(size(FTmat,1)/2):1:ceil(size(FTmat,1)/2)-1]);
    Rmat=((nx.^2)+(ny.^2)).^0.5;
    Angmat=atan2d(ny,nx)+0.01;
    corFT=FTmat*0;
    corFT=corFT+((Rmat<startFR).*imgFT);
    maxcount=min(floor(size(FTmat,2)/2),floor(size(FTmat,1)/2));
    maxcount=floor(maxcount*1.4*1);
    
    %loop over radii (spatial frequencies)
    for icount= startFR:maxcount-0
        %normalise moments
        P20C=CP20(icount)/CP00(icount);
        P22C=CP22(icount)/CP00(icount);

        % inverse profile and add ring to inverse operator
        corRing=((Rmat>=icount) & (Rmat<(icount+1))).*imgFT;
        corRing=corRing./(1+(P20C*cosd(2*Angmat))+(P22C*sind(2*Angmat)));
        corFT=corFT+corRing;
        
        %progress indication
        if ((rem(icount*100/maxcount,10) < (100/maxcount)) && (progFlag==true))
            p=round(icount*100/maxcount);
            message=sprintf(formspec,p);
            disp([backmsg message]);
        end
    end
    
    %Apply operator to image and convert back to real space
    corFT(ceil(size(corFT,1)/2)+1,ceil(size(corFT,2)/2)+1)=imgFT(ceil(size(corFT,1)/2)+1,ceil(size(corFT,2)/2)+1);
    corImg=abs(ifft2(fftshift(corFT)));

    %progress
    if (progFlag==true);
        disp(backmsg)
        fprintf('\b')
    end
end