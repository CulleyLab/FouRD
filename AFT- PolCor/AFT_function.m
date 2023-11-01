%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code and documentation can be found at
% https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform
%
% This routine uses a vector field of alignment directions using small
% sub-windows in the real space image to calculate an alignment order
% parameter.  Images should be grayscale images.
% Angles determined are oriented as follows:
%
%                            ^  180°
%                            |
%                            |
%                            ------>  90°
%                            |
%                            |
%                            v  0°
%
% Order Parameter value:   0 = Completely random alignment
%                          1 = Perfectly aligned
%
% All rights and permissions belong to:
% Patrick Oakes
% poakes@gmail.com
% 05/15/2015
%
% Citation:
% Cetera et al. Nat Commun 2014; 5:1-12
% http://www.ncbi.nlm.nih.gov/pubmed/25413675
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [av_ordermat,anglemat,exccentricitymat] = AFT_function(file, directory, parameters)

% load image
im = imread(fullfile(directory, file));
im = im2double(im);

% calculate angle vector field
[anglemat,exccentricitymat,pc,pr,vc2,ur2,excc] = AFT_anglemat(im, parameters);
% replace any NaNs in eccentricity with zeros
miscol=excc*0;
if (sum(isnan(excc)) > 0)
    miscol=or(miscol,isnan(excc));
    excc=fillmissing(excc,'spline');
end
% cap excentricity values to 1
if (max(excc)>1)
    miscol=or(miscol,double(excc>1));
    excc=max(excc.*(1-(excc>1)),(excc>1));
end
miscol=double(miscol);

% plots 
if parameters.figures == 1
    
    % vector field
    if (parameters.parproc == 0)
        figure
    else
        figure('Visible','off')
    end
    imshow(im, [])
    hold on
    if exist('pc','var') ~= 0
        for k = 1:numel(pc)    
            if (parameters.colvectors==1)
                quiver(pc(k),pr(k),vc2(k),ur2(k),0,'showarrowhead','off','linewidth',2,'Color',[excc(k),miscol(k),1-excc(k)])
            else
                quiver(pc(k),pr(k),vc2(k),ur2(k),0,'showarrowhead','off','linewidth',2,'Color',[1,1,0])
            end
        end
    end
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'], ['vectors_' file(1:end-4) '.tif']));
    close
    
    % angle heat map
    if (parameters.parproc == 0)
        figure('Position',[100 100 size(im,2)/3 size(im,1)/3]);
    else
        figure('Position',[100 100 size(im,2)/3 size(im,1)/3],'Visible','off');
    end
    imagesc(rad2deg(anglemat));
    hsv_nan = [[0,0,0];colormap('hsv')];
    set(gca,'visible','off')
    caxis([0,180]);
    colormap(hsv_nan);
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);
    %set black for NaNs only
    imagesc(max(rad2deg(anglemat),-180*isnan(anglemat)));
    caxis([-180/(size(hsv_nan,1)-1),180]);
    
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'], ['angle_heatmap_' file(1:end-4) '.tif']));
    close

    % eccentricity heat map
    if (parameters.parproc == 0)
        figure('Position',[200 100 size(im,2)/3 size(im,1)/3]);
    else
        figure('Position',[200 100 size(im,2)/3 size(im,1)/3],'Visible','off');
    end
    imagesc(exccentricitymat);
    %hsv_nan = [[0,0,0];colormap('hsv')];
    set(gca,'visible','off')
    caxis([0,1]);
    colormap(parula);
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);

    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'], ['excc_heatmap_' file(1:end-4) '.tif']));
    close
end

% calculate order parameter
av_ordermat = AFT_ordermat(anglemat, parameters);

end