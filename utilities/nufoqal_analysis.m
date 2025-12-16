function [error] = nufoqal_analysis(analysisFolder, resultsFolder, saveImg, params, onylOne, filename, nucleitag)
% Main function of NuFoQal analysis
% Read and associate the files
% Prepare the different images - denoise
% Drive the segmentations and quantifications
%
% Micka�l Garnier - UMR3664 - PICT@PASTEUR - 2025

%% Initialisation
% Get parameters and prepare folders
if iscell(analysisFolder)
    Loc_data = analysisFolder{1};
    Loc_data2 = analysisFolder{2};
    if onylOne
        %filename2 = filename{2};
        filename = filename{1};
    end
    two_sets = true;
else
    Loc_data = analysisFolder;
    two_sets = false;
end
Loc_resu = resultsFolder;

resXY = str2double(params{1});
resZ =  str2double(params{2});
nucleimaxradius = str2double(params{3});
spotsize =  eval(params{4});
toSkip = str2double(params{5});
Sig_info = params{6};
Sig_qual = eval(params{7});
use_cellpose = eval(params{8});
if isempty(params{9})
    minSigTh = nan(1);
else
    minSigTh = eval(params{9});
end
rad_dilate = eval(params{10});
if nucleitag
    nuclei_channel = params{12};
end

if two_sets
    fprintf('Data locations: ');
    cprintf('hyper', '%s\t%s', Loc_data, Loc_data2);
    if nucleitag
        fprintf('\n\t Nuclei images tag: %s', nuclei_channel);
    end
    fprintf('\nResults location: ');
    cprintf('hyper', '%s\n\n', Loc_resu);
else
    fprintf('Data locations: ');
    cprintf('hyper', '%s', Loc_data);
    if nucleitag
        fprintf('\n\t Nuclei images tag: %s', nuclei_channel);
    end
    fprintf('\nResults location: ');
    cprintf('hyper', '%s\n\n', Loc_resu);
end
fprintf('\t---\n\tParameters:\n');
fprintf('\tPixel size: %.3fx%.3fx%.3f\n', resXY, resXY, resZ);
fprintf('\tMaximal nuclei radius: %d', nucleimaxradius)
fprintf('\tSpot sizes:');
for ii = 1:length(spotsize)
    fprintf(' %.2f', spotsize(ii));
end
fprintf('\n\tMinimum variations (in stds): %.2f\n', minSigTh);
fprintf('\tSignal tag: %s\n', Sig_info);
fprintf('\tLow signal?: %d\n\t---\n', Sig_qual);
fprintf('\tAdd pixels around nuclei masks: %d\n\t---\n', rad_dilate);


if(~isfolder(Loc_resu))
	mkdir(Loc_resu);
end

% Set inter parameters and get files
error = 0;

cd(Loc_data);

files = dir;
files = files(3:end);

nbFiles = length(files);
tokill = zeros(nbFiles, 1);
for i_f = nbFiles:-1:1
    if files(i_f).isdir
        tokill(i_f) = 1;
    else
        try
            imfinfo(files(i_f).name);
        catch err
            fprintf('File %s: skipped !\n%s\n', files(i_f).name, err.message);
            tokill(i_f) = 1;
        end
    end
end
files(tokill == 1) = [];
nbFiles = length(files);

if two_sets
    cd(Loc_data2);
    files2 = dir;
    files2 = files2(3:end);
    
    nbFiles2 = length(files2);
    tokill = zeros(nbFiles2, 1);
    for i_f = nbFiles2:-1:1
        if files2(i_f).isdir
            tokill(i_f) = 1;
        else
            try
                imfinfo(files2(i_f).name);
            catch err
                fprintf('File %s: skipped !\n%s\n', files2(i_f).name, err.message);
                tokill(i_f) = 1;
            end
        end
    end
    files2(tokill == 1) = [];
    nbFiles2 = length(files2);
    
    if nbFiles ~= nbFiles2
        warning('Warnings:nbFilesnbFiles2', 'Numbers of files in both directory are not the same! Images with missing correspondance will be skipped!');
    end
end

% % Find files delimiter and basename positions
% fnames = struct2cell(files);
% fnames = fnames(1,:)';
% 
% delimiters = {'-', '_', ' '};
% deltest = zeros(nbFiles, 3);
% for id = 1:3
%     deltest(:, id) = cellfun(@(x) size(split(x, delimiters{id}), 1), fnames, 'UniformOutput', true);
% end
% [~, idel] = max(sum(deltest, 1));
% numparts = min(deltest(:, idel));

%% Locate signal tag images (and nuclei tag is required)
gfpimgs = zeros(nbFiles, 1);
if nucleitag
    nucsimgs = zeros(nbFiles ,1);
end
for i_f = 1:nbFiles
    % Check for signal tag
    if strfind(files(i_f).name, Sig_info)
        gfpimgs(i_f) = 1;
    elseif strfind(files(i_f).name, upper(Sig_info))
        gfpimgs(i_f) = 1;
    end

    % Check for nuclei tag
    if nucleitag
        if strfind(files(i_f).name, nuclei_channel)
            nucsimgs(i_f) = 1;
        elseif strfind(files(i_f).name, upper(nuclei_channel))
            nucsimgs(i_f) = 1;
        end
    end
end
nbGFP = sum(gfpimgs);
gfpimgs = find(gfpimgs);

if nucleitag
    nbNUCS = sum(nucsimgs);
    nucsimgs = find(nucsimgs);
    if nbNUCS ~= nbGFP
        warning('Warnings:NucsGfpNb', 'Numbers of signal tag and nuclei tag images are not the same! Images with missing correspondance will be skipped!');
    end
end
if two_sets
    if nucleitag
        decnucsimgs = zeros(nbFiles ,1);
    end
    decimgs = zeros(nbFiles2, 1);
    for i_f = 1:nbFiles2
        % Check for dec
        if strfind(files2(i_f).name, Sig_info)
            decimgs(i_f) = 1;
        elseif strfind(files2(i_f).name, upper(Sig_info))
            decimgs(i_f) = 1;
        end

        % Check for nuclei tag
        if nucleitag
            if strfind(files2(i_f).name, nuclei_channel)
                decnucsimgs(i_f) = 1;
            elseif strfind(files2(i_f).name, upper(nuclei_channel))
                decnucsimgs(i_f) = 1;
            end
        end
    end
    nbDec = sum(decimgs);
    % nbDec = sum(decimgs);
    if nbDec ~= nbGFP
        warning('Warnings:DecGfpNb', 'Numbers of quantification and segmentation signal tag images are not the same! Images with missing correspondance will be skipped!');
    end
    decimgs = find(decimgs);

    if nucleitag
        nbDECNUCS = sum(decnucsimgs);
        decnucsimgs = find(decnucsimgs);
        if nbDECNUCS ~= nbDec
            warning('Warnings:DecNucsDecNb', 'Numbers of deconvolved signal tag and deconvolved nuclei tag images are not the same! Images with missing correspondance will be skipped!');
        end
    end
end

%% Analysis loop on every files
for i_file = 1:nbGFP
    tic;
    cd(Loc_data);
    i_gfp = gfpimgs(i_file);
    if strcmp(files(i_gfp).name(end-3:end), '.stk')
        continue
    end
    
    if onylOne
        if ~strcmp(files(i_gfp).name, filename)
            continue
        end
    end

    infopos = strfind(files(i_gfp).name, Sig_info);
    if isempty(infopos)
        infopos = strfind(files(i_gfp).name, upper(Sig_info));
    end
    basename = files(i_gfp).name(1:infopos(1)-2);
    basename = erase(basename, 'corrected');

    if nucleitag
        % Get nuclei image
        for i_file2 = 1:nbNUCS
            if strfind(files(nucsimgs(i_file2)).name, basename)
                i_nucs = nucsimgs(i_file2);
                break;
            end
        end
    end
    
    if two_sets
        % Get dec image
        for i_file2 = 1:nbDec
            if strfind(files2(decimgs(i_file2)).name, basename)
                i_dec = decimgs(i_file2);
                break;
            end
        end
    end

    if two_sets & nucleitag
        % Get nuclei image
        for i_file2 = 1:nbNUCS
            if strfind(files2(decnucsimgs(i_file2)).name, basename)
                i_decnucs = decnucsimgs(i_file2);
                break;
            end
        end
    end

    fprintf('\nAnalysing files:');
    cprintf('Text*', '\t%s\t: %s', Sig_info, files(i_gfp).name);
    if two_sets
        cprintf('Text*', '\n\tSpots segmentation: %s', files2(i_dec).name);
    end
    if nucleitag
        cprintf('Text*', '\n\tNuclei segmentation: %s', files(i_nucs).name);
    end
    if two_sets & nucleitag
        cprintf('Text*', '\n\tDec nuclei segmentation: %s', files2(i_decnucs).name);
    end
    
    fprintf('\n\t> Reading data information and preparing output files. \n');
    %% File reading
    info = imfinfo(files(i_gfp).name);
    
    dims = uint16([info(1).Height info(1).Width length(info)]);
    
    switch info(1).BitDepth
        case 16
            type = 'uint16';
        case 8
            type = 'uint8';
        otherwise
            type = 'double';
    end
    
%     if ~noTrans
%         Trans = imread(files(i_trans).name);
%     end
    dims(3) = dims(3) - toSkip;
    GFP = zeros(dims(1), dims(2), dims(3), type);
    if two_sets
        DEC = zeros(dims(1), dims(2), dims(3), type);
    end
    if nucleitag
        NUCS = zeros(dims(1), dims(2), dims(3), type);
    end
    if two_sets & nucleitag
        DECNUCS = zeros(dims(1), dims(2), dims(3), type);
    end
    for iz = (1+toSkip):dims(3)
        GFP(:,:,iz) = imread(files(i_gfp).name, iz);
        if two_sets
            DEC(:,:,iz) = imread([Loc_data2, filesep, files2(i_dec).name], iz);
        end
        if nucleitag
            NUCS(:,:,iz) = imread(files(i_nucs).name, iz);
        end
        if two_sets & nucleitag
            DECNUCS(:,:,iz) = imread([Loc_data2, filesep, files2(i_decnucs).name], iz);
        end
    end

    %% STEP 0: Correct for background variations and denoise
    transformed = double(GFP);
    bckg_med = zeros(dims(3), 1);
    for i_z = 1:dims(3)
        slice = transformed(:,:,i_z);
        bckg_med(i_z) = median(slice(:));
    end

    for i_z = 1:dims(3)
        transformed(:,:,i_z) = transformed(:,:,i_z)./max(1, bckg_med(i_z));
    end

    transformed = Anscombe_forward(transformed);
    wienersize = max(min(spotsize), 3)*2;
    for i_z = 1:dims(3)
%         transformed(:,:,i_z) = medfilt2(transformed(:,:,i_z), [wienersize wienersize]);
        transformed(:,:,i_z) = wiener2(transformed(:,:,i_z), [wienersize wienersize]);
    end
    GFPmod = Anscombe_inverse_exact_unbiased(transformed);

    if nucleitag
        transformed = double(NUCS);
        bckg_med = zeros(dims(3), 1);
        for i_z = 1:dims(3)
            slice = transformed(:,:,i_z);
            bckg_med(i_z) = median(slice(:));
        end

        for i_z = 1:dims(3)
            transformed(:,:,i_z) = transformed(:,:,i_z)./max(1, bckg_med(i_z));
        end

        transformed = Anscombe_forward(transformed);
        wienersize = max(min(spotsize), 3)*2;
        for i_z = 1:dims(3)
    %         transformed(:,:,i_z) = medfilt2(transformed(:,:,i_z), [wienersize wienersize]);
            transformed(:,:,i_z) = wiener2(transformed(:,:,i_z), [wienersize wienersize]);
        end
        NUCSmod = Anscombe_inverse_exact_unbiased(transformed);
    end
    
    if two_sets
        transformed = double(DEC);
        bckg_med = zeros(dims(3), 1);
        for i_z = 1:dims(3)
            slice = transformed(:,:,i_z);
            bckg_med(i_z) = median(slice(:));
        end

        for i_z = 1:dims(3)
            transformed(:,:,i_z) = transformed(:,:,i_z)./max(1, bckg_med(i_z));
        end

        transformed = Anscombe_forward(transformed);
        wienersize = max(min(spotsize), 3)*2;
        for i_z = 1:dims(3)
    %         transformed(:,:,i_z) = medfilt2(transformed(:,:,i_z), [wienersize wienersize]);
            transformed(:,:,i_z) = wiener2(transformed(:,:,i_z), [wienersize wienersize]);
        end
        DECmod = Anscombe_inverse_exact_unbiased(transformed);
    end

    if two_sets & nucleitag
        transformed = double(DECNUCS);
        bckg_med = zeros(dims(3), 1);
        for i_z = 1:dims(3)
            slice = transformed(:,:,i_z);
            bckg_med(i_z) = median(slice(:));
        end

        for i_z = 1:dims(3)
            transformed(:,:,i_z) = transformed(:,:,i_z)./max(1, bckg_med(i_z));
        end

        transformed = Anscombe_forward(transformed);
        wienersize = max(min(spotsize), 3)*2;
        for i_z = 1:dims(3)
    %         transformed(:,:,i_z) = medfilt2(transformed(:,:,i_z), [wienersize wienersize]);
            transformed(:,:,i_z) = wiener2(transformed(:,:,i_z), [wienersize wienersize]);
        end
        DECNUCSmod = Anscombe_inverse_exact_unbiased(transformed);
    end
    
    %% STEP 1: NUCLEI SEGMENTATION (in 3D)
    if nucleitag
        if two_sets
            map3D = (mat2gray(DECNUCSmod) + mat2gray(GFPmod))./2;
        else
            map3D = double(NUCSmod);
        end
    else
        map3D = double(GFPmod);
    end

    try
        seg_gfp = main_nuclei_segmentation(map3D, use_cellpose, dims, nucleimaxradius, max(spotsize)+1, Sig_qual, rad_dilate, params(13:15));
    catch err
        cprintf('Error', err.message);
        continue
    end

    clear map3D

    if max(seg_gfp(:)) == 0
        cprintf('Error', 'No Nuclei have been found\n');
        outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.nucs');
        fid = fopen(outputfile, 'w+');
        fprintf(fid, 'No Nuclei were found');
        fclose(fid);

        continue
    end
    
    %% STEP 2: SPOTS SEGMENTATION (in 3D)
    val_name = strcat(Loc_resu, '/',files(i_file).name(1:end-4), '_seg_valid.tif');
    try
        if two_sets
            [spots, seg_gfp, seg_spots] = spotsSegmentation_checkwithraw_stats(GFPmod, DECmod, GFP, DEC, seg_gfp, dims, spotsize, resXY, resZ, 0.7, val_name, Sig_qual, minSigTh);
        else
            [spots, seg_gfp, seg_spots] = spotsSegmentation_stats(GFPmod, GFPmod, GFP, GFP, seg_gfp, dims, spotsize, resXY, resZ, 0.7, val_name, Sig_qual, minSigTh);
        end
    catch err
        cprintf('Error', err.message);
        continue
    end
    if size(spots,1) == 0
        cprintf('Error', 'No foci have been found\n');
        outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.nucs');
        fid = fopen(outputfile, 'w+');
        fprintf(fid, 'No foci have been found');
        fclose(fid);

        continue
    end
    

    nbSpots = size(spots,1);
    if saveImg
        seg_name = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '_segSpots.tif');
        imwrite(uint16(seg_spots(:,:,1)), seg_name, 'writemode', 'overwrite');
        for i = 2:dims(3)
            imwrite(uint16(seg_spots(:,:,i)), seg_name, 'writemode', 'append');
        end
        seg_name = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '_segNuc.tif');
        imwrite(uint16(seg_gfp(:,:,1)), seg_name, 'writemode', 'overwrite');
        for i = 2:dims(3)
            imwrite(uint16(seg_gfp(:,:,i)), seg_name, 'writemode', 'append');
        end
    end

    %% STEP 3: NUCLEI STATISTICS (in 3D)
    try
        [nuclei] = measureInfosNuclei(GFP, GFPmod, seg_gfp, seg_spots, spots, saveImg, strcat(Loc_resu, '/',files(i_file).name(1:end-4), '_segSpotsSpheres.tif'));
    catch err
        cprintf('Error', err.message);
        continue
    end

    %% WRITE OUTPUT  NUCLEI % --------------------------------------------------------- %
	% Nucleus: 1> Intensity					Output: 1> Nucleus Volume                   %
	% Nucleus: 2> Intensity background				2> Nucleus intensity                %
	% Nucleus: 3> Number of foci					3> Nucleus intensity background     %
	% Nucleus: 4> Volume							4> Nucleus number of spots      	%
	% Nucleus: 5> X         						5> ToKeep                           %
    % Nucleus: 6> Y                                 6> X    7> Y    8> Z                %
    % Nucleus: 7> Z                                                                     %
    % --------------------------------------------------------------------------------- %
	
    outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.nucs');
    fid = fopen(outputfile, 'w+');
    fprintf(fid, 'Nucleus intensity\tNucleus background\tNucleus Volume\tNucleus spot number\tTo keep\tX\tY\tZ\n');

    bckg_int = median(double(GFP(:)));
	nbNucs = size(nuclei,1);
    
    tokeepnuc = nuclei(:,4) > 999; % Hard-coded to avoir smaller chunk of signal (dust and whatnot)
    
    for i = 1:nbNucs
		%              I Ibg   V  nbSpots  keep
        fprintf(fid, '%d\t%d\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f', ...
            nuclei(i, 1)-bckg_int*nuclei(i, 4), nuclei(i, 2)-bckg_int*nuclei(i, 4), nuclei(i, 4), nuclei(i, 3), tokeepnuc(i), nuclei(i,5), nuclei(i,6), nuclei(i,7));
        if i < nbNucs
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    %% WRITE OUTPUT  SPOTS % --------------------------------------------------------------------------------------------------------------------------------- %
	% Spots: 1> X						Nucleus: 1> Intensity					Output: 1> X spot							8> Relative distance to border		%
	%		 2> Y						Nucleus: 2> Intensity background				2> Y spot							9> Nucleus ID						%
	%		 3> Z						Nucleus: 3> Number of foci						3> Z spot							A> Nucleus intensity				%
	%		 4> Volume					Nucleus: 4> Volume								4> Volume spot						B> Nucleus intensity background		%
	%		 5> Intensity                                   							5> Intensity spot					C> Nucleus Volume					%
	%		 6> Nucleus ID																6> Relative spot intensity			D> Nucleus number of spots			%
	%		 7> distance to border		8> Sphericity (if low)							7> Distance to nucleus border		E> ToKeep		F> tElongation		%
	% --------------------------------------------------------------------------------------------------------------------------------------------------------- %
	
    outputfile = strcat(Loc_resu, '/',files(i_gfp).name(1:end-4), '.spots');
    fid = fopen(outputfile, 'w+');
    fprintf(fid, 'X\tY\tZ\tVolume\tIntensity\tRelative intensity\tDistance to border\tRelative distance to border\tNucleus ID\tNucleus intensity\tNucleus background\tNucleus Volume\tNucleus spot number\tTo keep\tElongation\n');

    bckg_int = median(double(GFP(:)));
	nuclei_radii = (3/(4*pi)*nuclei(:,4)).^(1/3);
    
    for i = 1:nbSpots
        relInt = (spots(i,5) -bckg_int*spots(i,4))/(nuclei(spots(i,6), 1) - nuclei(spots(i,6), 4) * bckg_int);
		%              X     Y     Z     V      I     Ir    D     Dr   nuc  I Ibg   V  nbSpots  keep   sphericity
        fprintf(fid, '%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.4f\t%.4f\t%.2f\t%d\t%d\t%d\t%.0f\t%.0f\t%.0f\t%d', ...
            spots(i,1), spots(i,2), spots(i,3), spots(i,4), spots(i,5)-spots(i,4)*bckg_int, relInt, spots(i,7), spots(i,7)/nuclei_radii(spots(i,6)),...
            spots(i,6), nuclei(spots(i,6), 1)-bckg_int*nuclei(spots(i,6), 4), nuclei(spots(i,6), 2)-bckg_int*nuclei(spots(i,6), 4), nuclei(spots(i,6), 4), nuclei(spots(i,6), 3), tokeepnuc(spots(i,6)), spots(i,8));
        if i < nbSpots
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    toc;
end
cprintf([1 0.5216 0], '\nqFOCI analysis finished !\n-----------------------------\n');

