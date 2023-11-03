% This code is a modified version of the 'Alignment by Fourier Transform' package associated with the paper
% "A Workflow for Rapid Unbiased Quantification of Fibrillar Feature Alignment in Biological Images" Front. Comput. Sci 14 October 2021, Stefania Marcotti et.al.
% DOI: https://doi.org/10.3389/fcomp.2021.745831
% Original software avalible at https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform
% All rights and permissions belong to: Patrick Oakes poakes@gmail.com
% The modifications support intigration into the FouRD package submitied with the paper
% "Characterisation and correction of polarisation effects in fluorescently labelled fibres", Nandini Aggarwi et.al.
% submitted to the Journal of Microscopy. 3/11/2023

function av_ordermat = AFT_ordermat(anglemat, parameters)

% parameters and initialisation
st = parameters.st;
ordermat = NaN(size(anglemat));

% for each neighbourhood of size 2st+1 vectors
for i = st+1:size(anglemat,1)-st
    for j = st+1:size(anglemat,2)-st
        
        % compare ref vector to neighbourhood
        temp = anglemat(i-st:i+st,j-st:j+st);
        temp2 = repmat(anglemat(i,j),2*st+1,2*st+1);
        comp = cos(temp-temp2);
        comp = comp.*comp;
        
        % remove central value (comparison ref with itself)
        idx_centre = round(size(comp,1)/2);
        comp(idx_centre, idx_centre) = NaN;
        
        % calculate order parameter
        ordermat(i,j) = 2*(mean(comp(:))-.5);
        
        clear temp temp2 comp
    end
end

% calculate median order parameter for image
av_ordermat = median(ordermat(:));

end