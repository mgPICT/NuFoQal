function [Sout] = illumination_correction(img_in, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    dims = size(img_in);

    options = struct;
    %% options that may be specified by the user
    options.lambda_vreg             = 6;   % default value = 6
    options.lambda_zero             = 0.5;   % default value = 0.5
    options.max_lbfgs_iterations    = 500;   % default value = 500
    options.q_percent               = 0.25;   % default value = 0.25
    options.image_size              = dims(1:2);
    options.folder_source           = [];
    options.folder_destination      = [];
    options.filenames               = {};
    options.num_images_provided     = [];
    options.bit_depth               = 2^16;   % specified as maximum integer: 2^8, 2^12, 2^16
    options.correction_mode         = 0;   % 0 ='illumination preserving' (default), 1='zero-light_perserving', or 2='direct'

    %% internal options, should not be reset by user without expert knowledge
    options.target_num_pixels     	= 9400;
    options.working_size            = [];
    options.number_of_quantiles     = 200;
    %options.lambda_barr             = [];
    

    [S, options] = cdr_loadImages(img_in, options);

    [model, options] = cdr_cidreModel(S, options);

    %% Apply correction from cdr_correct on images on the go (skip read/write part)
    Sout = zeros(dims);

    for z = 1:options.num_images_provided
        I = double(img_in(:,:,z));

        % check which type of correction we want to do
        switch options.correction_mode
            case 0  %'intensity range _preserving'
                Icorrected = ((I - model.z)./model.v) * mean(model.v(:))  + mean(model.z(:));

            case 1 % 'zero-light_preserving'
                Icorrected = ((I - model.z)./model.v) * mean(model.v(:));

            case 2 %'direct'    
                Icorrected = ((I - model.z)./model.v);

            otherwise
                error('CIDRE:correction', 'Unrecognized correction mode: %s', lower(options.correction_mode));
        end
        Sout(:,:,z) = Icorrected;
    end


end

