function labels = main_nuclei_segmentation(pixmap, use_cellpose, dims, maxradius, spotsizedef, Sig_qual, rad_dilate, cp_params)

    %% Grab parameters and prepare input/output
    fprintf('3D segmentation in progress\n');
    spotsize = round(sqrt(2)*spotsizedef);
    
    output = double(pixmap);
    sumproj = sum(output, 3);

    if Sig_qual
        bckgproj = sumproj;
    else
        bckgproj = imtophat(sumproj, strel('disk', spotsize*4, 0));
    end

    %% 2D localisation of nuclei
    if use_cellpose
        labels = nuclei_2D_localisation_cellpose(bckgproj, cp_params{1}, eval(cp_params{2}), eval(cp_params{3}));
        labels = bwlabel(labels > 0, 4);
    else
        labels = nuclei_2D_localisation_legacy(bckgproj, dims, spotsize);
    end

    nbObject = max(labels(:));
    fprintf('2D localisations: %d nuclei found\n', nbObject);

    %% 3D nucleus-wise segmentation
    labels = nuclei_3D_segmentation(pixmap, labels, dims, nbObject, maxradius, rad_dilate);

end