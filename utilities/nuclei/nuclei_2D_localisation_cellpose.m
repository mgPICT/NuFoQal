function labels = nuclei_2D_localisation_cellpose(map2D, cp_model, cellDiam, useEns)
%% nuclei_2D_localisation_cellpose(map2D, cellDiam, useEns)
%
% Locate and segment nuclei in 2D from map2D image using Cellpose
% Remove objects on edges and small chunks of pixels
% Relabel everything after cleaning
%
% Inputs:
%   map2D - image to segment
%   cellDiam - approximate diameter of nuclei to detect in pixels
%   useEns - boolean, use all versions of trained cyto models of cellpose

    cp = cellpose('Model', cp_model, UseEnsemble=useEns);
    labels = segmentCells2D(cp, map2D, ImageCellDiameter=cellDiam);
    
    labels = imclearborder(labels, 8);
    conts = imdilate(edge(labels, 'Sobel', 0), ones(3));
    labels(conts) = 0;
    labels = bwlabel(labels > 0, 4);
end
