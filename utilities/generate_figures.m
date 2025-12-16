function [errorout] = generate_figures(resultsFolder, params)
%% generate_figures
%  Read qFOCI2 .nucs files and generate nice figures from the results. 
% 
%% Output
%  Figures designed with Antoine Hocher on R displaying:
%      > the histogram of the number of foci per cell
%      > the pdfs of the percentage of intensity contained in the foci
%      > the pdfs of the total intensity contained in the foci
%      > the pdfs of the foci position in the cell (distance to border)
%  for each cell line (WT and mutants)
%
%  Bonus figures:
%      > Boxplots per yeastLines of intensities (absolutes and relatives)
%      > Pdfs of relatives and absolutes intensities of spots
%
% Micka�l Garnier - UMR3664 - PICT@PASTEUR - 2016

%% Initialisation
% Get parameters and prepare folders
Loc_resu = resultsFolder;
resXY = str2double(params{1});
resZ =  str2double(params{2});
Sig_info = params{6};
Yeast_info = params{11};
if(~isfolder(Loc_resu))
	errorout = 'The selected result folder does not exist, you need to run the analysis first !';
    return
end

% Set inter parameters and get files
errorout = 0;

cd(Loc_resu);

% Manually select the files to analyse
[files, nucsPath] = uigetfile([Loc_resu, filesep, '*.spots'], 'Select the .spots files to check', 'MultiSelect', 'on');

% Enter path to figures
graphDir = uigetdir(Loc_resu, 'Select the folder to save the figures!');

nbFiles = size(files, 2);
if ~iscell(files)
    files = {files};
    nbFiles = 1;
end

yeastTypes = cell([1, 128]);
fileYeast = zeros(nbFiles, 1);

cd(nucsPath);

i_yeast = 0;
for i_f = 1:nbFiles
    id_type_beg = strfind(files{i_f}, Yeast_info);
    id_type_end = id_type_beg+strfind(files{i_f}(id_type_beg:end), '_');
    id_type_end = id_type_end(1);
    id_type_end2 = id_type_beg+strfind(files{i_f}(id_type_beg:end), '-');
    id_type_end2 = id_type_end2(1);
    id_type_end = min([id_type_end, id_type_end2]);
    offset = 2;
    if isempty(id_type_end)
        id_type_end = strfind(files{i_f}, upper(Sig_info));
        offset = 4;
    end
    if isempty(id_type_end)
        id_type_end = strfind(files{i_f}, Sig_info);
        offset = 4;
    end
    id_type_beg = id_type_beg(1);
    id_type_end = id_type_end(end);
    
    yeastID = files{i_f}(id_type_beg:id_type_end-offset);
    
    if isempty(yeastTypes)
        yeastTypes{i_yeast +1} = yeastID;
        i_yeast = i_yeast +1;
        fileYeast(i_f) = i_yeast;
    else
        isnew = true;
        for i_y = 1:i_yeast
            if strcmp(yeastID, yeastTypes{i_y})
                isnew = false;
                fileYeast(i_f) = i_y;
                break
            end
        end
        if isnew
            yeastTypes{i_yeast +1} = yeastID;
            i_yeast = i_yeast +1;
            fileYeast(i_f) = i_yeast;
        end
    end
end
yeastTypes = yeastTypes(~cellfun(@isempty, yeastTypes));
[yeastTypes, ordyt] = sort(yeastTypes);
[~, rev_ordyt] = sort(ordyt);
fileYeast = rev_ordyt(fileYeast)';

disp(yeastTypes);

nbYeastTypes = size(yeastTypes, 2);

if nbYeastTypes == 1 % Ask user for specific conditions
    answer = questdlg('Only one yeast strain was found, do you want to define specific filters ?');
    if strcmp(answer, 'Yes')
        nbYeastTypes = str2double(inputdlg('Number of groups'));
        prompt = cell([nbYeastTypes, 1]);
        for iy = 1:nbYeastTypes
            prompt{iy} = sprintf('Group %d', iy);
        end
        yeastTypes = inputdlg(prompt)';
        
        % Regroup files
        for i_f = 1:nbFiles
            fileYeast(i_f) = 0;
            for iy = 1:nbYeastTypes
                if ~isempty(strfind(files{i_f}, yeastTypes{iy}))
                    fileYeast(i_f) = iy;
                    continue
                end
            end
        end
        
        
        yeastTypes = yeastTypes(~cellfun(@isempty, yeastTypes));
        [yeastTypes, ordyt] = sort(yeastTypes);
        [~, rev_ordyt] = sort(ordyt);
        fileYeast = rev_ordyt(fileYeast)';

        disp(yeastTypes);
    end
end

%% Read information and gather data from results files
nbFociHist = zeros(nbYeastTypes, 128);

cumulInt = cell([nbYeastTypes 3]);
cumulDists = cell([nbYeastTypes 5]); % Dist2edge / RelDist2center / Dist2center / Dist2edge_norad / RelDist2center_norad
cumulSpotsNb = cell([nbYeastTypes 1]);
cumulNucleusID = cell([nbYeastTypes 1]);

cumulVol = cell([nbYeastTypes 2]);

cumulNucFeats = cell([nbYeastTypes 3]);

emptynucs = zeros(nbYeastTypes, 1);

% Generate big measurement file
fid = fopen(fullfile(Loc_resu, 'all_measurements.txt'), 'w');
fprintf(fid, 'Image\tStrain\tX\tY\tZ\tVolume\tVolume (µm)\tIntensity\tRelative intensity\tDistance to border (µm)\tDistance to nucleus center (µm)\tRelative distance to center\tDistance to border (no radius) (µm)\tRelative distance to center (no radius) (µm)\tNucleus ID images\tNucleus ID unique\tNucleus intensity\tNucleus background\tNucleus volume\tNucleus volume (µm3)\tNucleus spot number\tTo keep\tElongation\n');

Nucleus_ID_offset = 0;
for i_file = 1:nbFiles
    imgName = files{i_file}(1:strfind(files{i_file}, '.')-1);
    strainID = yeastTypes{fileYeast(i_file)};
    [X, Y, Z, Volume, Intensity,Relativeintensity,Distancetoborder,        ...
        ~, Nucleus_ID, Nucleusintensity, Nucleusbackground, NucleusVolume, ...
        Nucleusspotnumber, ToKeep, Elongation] = importSpots(files{i_file});
    
    nucsfile = replace(files{i_file}, '.spots', '.nucs');
    [~,~,~,~,~, Xn, Yn, Zn] = importNucs(nucsfile);

    distance2center = sqrt((resXY*(X-Xn(Nucleus_ID))).^2 + (resXY*(Y-Yn(Nucleus_ID))).^2 + (resZ*(Z-Zn(Nucleus_ID))).^2);
    relativedistance2center = distance2center./(Distancetoborder+distance2center);

    distance2border_norad = Distancetoborder - (Volume*resXY*resXY*resZ.*3/4/pi).^(1/3);
    distance2border_norad(distance2border_norad < 0) = 0;
    distance2center_norad = distance2center + (Volume*resXY*resXY*resZ.*3/4/pi).^(1/3);
    relativedistance2center_norad = distance2center_norad./(distance2center_norad+distance2border_norad);

    nbSpots = size(X, 1);
    for ispot = 1:nbSpots
        %             img yat  X   Y  Z   V    Vum  I    Ir     d2b   d2c d2cr d2bnr d2brnr Ni Ni+  nI  nBg  nV nVum nSn  tk  elong
        fprintf(fid, '%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.2f\t%.4f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f', ...
            imgName, strainID, X(ispot), Y(ispot), Z(ispot), Volume(ispot), Volume(ispot)*resXY^2*resZ, ...
            Intensity(ispot), Relativeintensity(ispot), Distancetoborder(ispot), distance2center(ispot), ...
            relativedistance2center(ispot), distance2border_norad(ispot), relativedistance2center_norad(ispot), ...
            Nucleus_ID(ispot), Nucleus_ID(ispot)+Nucleus_ID_offset, Nucleusintensity(ispot), Nucleusbackground(ispot), ...
            NucleusVolume(ispot), NucleusVolume(ispot)*resXY^2*resZ, Nucleusspotnumber(ispot), ToKeep(ispot), Elongation(ispot));
        if ispot < nbSpots
            fprintf(fid, '\n');
        end
    end
    
    Nucleus_ID = Nucleus_ID+Nucleus_ID_offset;
	
    if i_file < nbFiles
        fprintf(fid, '\n');
        Nucleus_ID_offset = Nucleus_ID_offset + max(Nucleus_ID);
    end
end
fclose(fid);

% Read and compile data
Nucleus_ID_offset = 0;
for i_file = 1:nbFiles
    strainID = fileYeast(i_file);
    
    % Read .spots file
    [X, Y, Z, Volume, Intensity,Relativeintensity,Distancetoborder,        ...
        ~, Nucleus_ID, Nucleusintensity, Nucleusbackground, NucleusVolume, ...
                       Nucleusspotnumber, ToKeep, ~] = importSpots(files{i_file});
                   
    % Read .nucs file
    nucsfile = replace(files{i_file}, '.spots', '.nucs');
    [Nucleusintensity_nucs,Nucleusbackgroundintensity_nucs,NucleusVolume_nucs,Nucleusspotnumber_nucs,Tokeep_nucs, Xn, Yn, Zn] = importNucs(nucsfile);
    
    % Check files integrity !
    for ii = 1:size(Nucleus_ID, 1)
        if Nucleusintensity(ii) == Nucleusintensity_nucs(Nucleus_ID(ii))
            continue
        else
            errorout = 'Missmatch between nucleus intensities in *.spots and *.nucs files !';
            return
        end
    end
    
    Nucleus_ID_raw = Nucleus_ID;
    Nucleus_ID = Nucleus_ID+Nucleus_ID_offset;
    Nucleus_ID_offset = Nucleus_ID_offset + max(Nucleus_ID);
                   
    X(~ToKeep) = [];
    Y(~ToKeep) = [];
    Z(~ToKeep) = [];
    Intensity(~ToKeep) = [];
    Relativeintensity(~ToKeep) = [];
    Distancetoborder(~ToKeep) = [];
    %Relativedistancetoborder(~ToKeep) = [];
    Nucleus_ID_raw(~ToKeep) = [];
    Nucleus_ID(~ToKeep) = [];
    Nucleusintensity(~ToKeep) = [];
    Nucleusbackground(~ToKeep) = [];
    Nucleusspotnumber(~ToKeep) = [];
    
    Volume(~ToKeep) = [];
    NucleusVolume(~ToKeep) = [];
                   
    if max(Relativeintensity)>1
        Relativeintensity = Intensity./Nucleusintensity;
    end
    
    Nucleusintensity_nucs(~Tokeep_nucs) = [];
    Nucleusbackgroundintensity_nucs(~Tokeep_nucs) = [];
    NucleusVolume_nucs(~Tokeep_nucs) = [];
    Nucleusspotnumber_nucs(~Tokeep_nucs) = [];
    
    % Store information
    spotsNumH = histcounts(Nucleusspotnumber, 1:128);
    spotsNumH = spotsNumH./(1:127);
    nbFociHist(strainID,1:127) = nbFociHist(strainID,1:127) + spotsNumH;
    
    Intensity(Intensity == 0) = 1;
    Relativeintensity(Relativeintensity == 0) = 0.0001;
    cumulInt{strainID, 1} = [cumulInt{strainID, 1}; Intensity];
    cumulInt{strainID, 2} = [cumulInt{strainID, 2}; Relativeintensity];
    cumulInt{strainID, 3} = [cumulInt{strainID, 3}; Nucleusbackground];
    
    cumulSpotsNb{strainID} = [cumulSpotsNb{strainID}; Nucleusspotnumber];
    cumulNucleusID{strainID} = [cumulNucleusID{strainID}; Nucleus_ID];
    
    % Distancetoborder(Distancetoborder == 0) = 1;
    % Relativedistancetoborder(Relativedistancetoborder == 0) = 0.0001;

    distance2center = sqrt((resXY*(X-Xn(Nucleus_ID_raw))).^2 + (resXY*(Y-Yn(Nucleus_ID_raw))).^2 + (resZ*(Z-Zn(Nucleus_ID_raw))).^2);
    relativedistance2center = distance2center./(Distancetoborder+distance2center);

    distance2border_norad = Distancetoborder - (Volume*resXY*resXY*resZ.*3/4/pi).^(1/3);
    distance2border_norad(distance2border_norad < 0) = 0;
    relativedistance2center_norad = distance2center./(distance2center+distance2border_norad);

    cumulDists{strainID, 1} = [cumulDists{strainID, 1}; Distancetoborder];
    cumulDists{strainID, 2} = [cumulDists{strainID, 2}; relativedistance2center];
    cumulDists{strainID, 3} = [cumulDists{strainID, 3}; distance2center];
    cumulDists{strainID, 4} = [cumulDists{strainID, 4}; distance2border_norad];
    cumulDists{strainID, 5} = [cumulDists{strainID, 5}; relativedistance2center_norad];
    
    cumulVol{strainID, 1} = [cumulVol{strainID, 1}; Volume];
    cumulVol{strainID, 2} = [cumulVol{strainID, 2}; NucleusVolume];
    
    cumulNucFeats{strainID, 1} = [cumulNucFeats{strainID, 1}; Nucleusintensity_nucs];
    cumulNucFeats{strainID, 2} = [cumulNucFeats{strainID, 2}; Nucleusbackgroundintensity_nucs];
    cumulNucFeats{strainID, 3} = [cumulNucFeats{strainID, 3}; NucleusVolume_nucs];
    
    % Add empty nuclei
    emptynucs(strainID) = emptynucs(strainID) + sum(Nucleusspotnumber_nucs == 0);
end

% Grab maximum number of Foci per lines and intensity maxima
nbFociStrain = zeros(nbYeastTypes, 1);
maxIntS = zeros(nbYeastTypes, 1);
maxIntR = zeros(nbYeastTypes, 1);
maxDist = zeros(nbYeastTypes, 1);
for i = 1:nbYeastTypes
    nbFociStrain(i) = find(nbFociHist(i,:), 1, 'last');
    maxIntS(i) = max(cumulInt{i, 1});
    maxIntR(i) = max(cumulInt{i, 2});
    maxDist(i) = max(cumulDists{i, 1});
end
maxIS = max(maxIntS);
maxIR = max(maxIntR);
maxDs = max(maxDist);

% Compute pdfs
nbFociMax = max(nbFociStrain);%min(6, round(median(nbFociStrain)));
nbFociHist = cat(2, nbFociHist(:, 1:nbFociMax+1), sum(nbFociHist(:, nbFociMax+2:end), 2));

xi_pdf_intS = 0:maxIS/99:maxIS;
xi_pdf_intR = 0:maxIR/99:maxIR;
xi_pdf_dist = 0:maxDs/99:maxDs;
xi_pdf_distR = 1/100:1/100:1;

pdfIntS = cell([nbYeastTypes, nbFociMax]);
pdfIntR = cell([nbYeastTypes, nbFociMax]);
pdfDist = cell([nbYeastTypes, nbFociMax]);
pdfDistR = cell([nbYeastTypes, nbFociMax]);

pdfIntSum = cell([nbYeastTypes, nbFociMax]);
pdfIntPercent = cell([nbYeastTypes, nbFociMax]);

for i_y = 1:nbYeastTypes
    for inb = 1:nbFociMax
        id_inb = cumulSpotsNb{i_y} == inb;
        if sum(id_inb) == 0
            pdfIntR{i_y, inb} = zeros(size(xi_pdf_intS));
            pdfIntS{i_y, inb} = zeros(size(xi_pdf_intR));
            pdfDist{i_y, inb} = zeros(size(xi_pdf_dist));
            pdfDistR{i_y, inb} = zeros(size(xi_pdf_distR));
            pdfIntSum{i_y, inb} = zeros(size(xi_pdf_intS));
            pdfIntPercent{i_y, inb} = zeros(size(xi_pdf_intR));
            continue
        end
        pdfIntS{i_y, inb} = ksdensity(cumulInt{i_y, 1}(id_inb), xi_pdf_intS, 'support', [-0.01 maxIS*1.01]);
        pdfIntS{i_y, inb} = pdfIntS{i_y, inb}./sum(pdfIntS{i_y, inb});
        pdfIntR{i_y, inb} = ksdensity(cumulInt{i_y, 2}(id_inb), xi_pdf_intR, 'support', [-0.01 maxIR*1.01]);
        pdfIntR{i_y, inb} = pdfIntR{i_y, inb}./sum(pdfIntR{i_y, inb});
        pdfDist{i_y, inb} = ksdensity(cumulDists{i_y, 1}(id_inb), xi_pdf_dist, 'support', [-0.01 maxDs*1.01]);
        pdfDist{i_y, inb} = pdfDist{i_y, inb}./sum(pdfDist{i_y, inb});
        pdfDistR{i_y, inb} = ksdensity(cumulDists{i_y, 2}(id_inb), xi_pdf_distR, 'support', [-0.01 1.01]);
        pdfDistR{i_y, inb} = pdfDistR{i_y, inb}./sum(pdfDistR{i_y, inb});
        
        idNuclei = cumulNucleusID{i_y}(id_inb);
        [b, ~, n] = unique(idNuclei);
        
        IntS = cumulInt{i_y, 1}(id_inb);
        IntR = cumulInt{i_y, 2}(id_inb);
        
        IntSumNucs = [];
        IntPercent = [];
        for i_n = 1:length(b)
            IntSumNucs = [IntSumNucs; IntS(n == i_n)]; %#ok<AGROW>
            IntPercent = [IntPercent; IntR(n == i_n)]; %#ok<AGROW>
        end
        
        pdfIntSum{i_y, inb} = ksdensity(IntSumNucs, xi_pdf_intS, 'support', 'positive');
        pdfIntSum{i_y, inb} = pdfIntSum{i_y, inb}./sum(pdfIntSum{i_y, inb});
        pdfIntPercent{i_y, inb} = ksdensity(IntPercent, xi_pdf_intR, 'support', 'positive');
        pdfIntPercent{i_y, inb} = pdfIntPercent{i_y, inb}./sum(pdfIntPercent{i_y, inb});
    end
end

% Add empty ncus to histograms!
for i_y = 1:nbYeastTypes
    nbFociHist(i_y,:) = [emptynucs(i_y),  nbFociHist(i_y,1:nbFociMax+1)];
end

%% Create nice figures
cmap_yeast = jet(nbYeastTypes)./2;
cmap_nb = colorcube(nbFociMax+8).*(2/3);
h_all = cell([nbYeastTypes, 5]);

% Plot distributions, grouped by yeast strains and colored by number of
% spots per nuclei
for i_y = 1:nbYeastTypes
    legendstr = cell([1 nbFociStrain(i_y)]);
    h_all{i_y, 1} = figure('Visible', 'off');
    h_all{i_y, 2} = figure('Visible', 'off');
    h_all{i_y, 3} = figure('Visible', 'off');
    h_all{i_y, 4} = figure('Visible', 'off');
    h_all{i_y, 5} = figure('Visible', 'off');
    for inb = nbFociStrain(i_y):-1:1
        figure(h_all{i_y,1}); hold on, plot(xi_pdf_intS, pdfIntS{i_y, inb}, 'Color', cmap_nb(inb,:), 'LineWidth', 2);
        figure(h_all{i_y,2}); hold on, plot(xi_pdf_intR, pdfIntR{i_y, inb}, 'Color', cmap_nb(inb,:), 'LineWidth', 2);
        figure(h_all{i_y,3}); hold on, plot(xi_pdf_dist, pdfDist{i_y, inb}, 'Color', cmap_nb(inb,:), 'LineWidth', 2);
        figure(h_all{i_y,4}); hold on, plot(xi_pdf_distR, pdfDistR{i_y, inb}, 'Color', cmap_nb(inb,:), 'LineWidth', 2);
        if inb > 1
            legendstr{inb} = sprintf('%d spots', inb);
        else
            legendstr{inb} = sprintf('%d spot', inb);
        end
    end
    figure(h_all{i_y,1}); xlabel('Intensity'); ylabel('Probability'); title(yeastTypes{i_y}, 'interpreter', 'none'); legend(flip(legendstr));
    figure(h_all{i_y,2}); xlabel('Relative Intensity'); ylabel('Probability'); title(yeastTypes{i_y}, 'interpreter', 'none'); legend(flip(legendstr));
    figure(h_all{i_y,3}); xlabel('Distance'); ylabel('Probability'); title(yeastTypes{i_y}, 'interpreter', 'none'); legend(flip(legendstr));
    figure(h_all{i_y,4}); xlabel('Relative Distance'); ylabel('Probability'); title(yeastTypes{i_y}, 'interpreter', 'none'); legend(flip(legendstr));
end

% Plot histogram of nb Spots per yeast strain
for i_y = 1:nbYeastTypes
    figure(h_all{i_y,5}); bar(0:nbFociStrain(i_y), nbFociHist(i_y,1:nbFociStrain(i_y)+1), 'FaceColor', cmap_yeast(i_y,:));
    xlabel('Nb spots'); ylabel('Count'); title(yeastTypes{i_y}, 'interpreter', 'none');
end

% Write figs
for i_y = 1:nbYeastTypes
    set(h_all{i_y,1}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
    saveas(h_all{i_y,1}, [graphDir, filesep, yeastTypes{i_y}, '_intensity.pdf']);
    close(h_all{i_y,1});
    
    set(h_all{i_y,2}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
    saveas(h_all{i_y,2}, [graphDir, filesep, yeastTypes{i_y}, '_relative_intensity.pdf']);
    close(h_all{i_y,2});
    
    set(h_all{i_y,3}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
    saveas(h_all{i_y,3}, [graphDir, filesep, yeastTypes{i_y}, '_distance.pdf']);
    close(h_all{i_y,3});
    
    set(h_all{i_y,4}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
    saveas(h_all{i_y,4}, [graphDir, filesep, yeastTypes{i_y}, '_relative_distance.pdf']);
    close(h_all{i_y,4});
    
    set(h_all{i_y,5}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
    saveas(h_all{i_y,5}, [graphDir, filesep, yeastTypes{i_y}, '_histogram.pdf']);
    close(h_all{i_y,5});
end


% Generate boxplots
boxIntSum = [];
boxIntRel = [];
% boxIntBg = [];
boxVol = [];
boxVolRel = [];
boxD2B = [];
boxD2CRel = [];
boxD2C = [];
boxD2B_norad = [];
boxD2CRel_norad = [];
boxNbFoci = [];
groups = cell([0 0]);

boxIntNuc = [];
boxBgNuc = [];
boxVolNuc = [];
groupsnucs = cell([0 0]);
for i = 1:nbYeastTypes
    boxIntSum = [boxIntSum; cumulInt{i,1}]; %#ok<AGROW>
    boxIntRel = [boxIntRel; cumulInt{i,2}]; %#ok<AGROW>
%     boxIntBg = [boxIntBg; cumulInt{i,3}];
    boxVol = [boxVol; cumulVol{i,1}]; %#ok<AGROW>
    boxVolRel = [boxVolRel; cumulVol{i,1}./cumulVol{i,2}]; %#ok<AGROW>

    boxD2B = [boxD2B; cumulDists{i, 1}]; %#ok<AGROW>
    boxD2CRel = [boxD2CRel; cumulDists{i, 2}]; %#ok<AGROW>
    boxD2C = [boxD2C; cumulDists{i, 3}]; %#ok<AGROW>
    boxD2B_norad = [boxD2B_norad; cumulDists{i, 4}]; %#ok<AGROW>
    boxD2CRel_norad = [boxD2CRel_norad; cumulDists{i, 5}]; %#ok<AGROW>

    boxNbFoci = [boxNbFoci; cumulSpotsNb{i}]; %#ok<AGROW>

    groups = [groups; repmat(yeastTypes(i), size(cumulInt{i,1}))]; %#ok<AGROW>
    
    boxIntNuc = [boxIntNuc; cumulNucFeats{i,1}]; %#ok<AGROW>
    boxBgNuc = [boxBgNuc; cumulNucFeats{i,2}]; %#ok<AGROW>
    boxVolNuc = [boxVolNuc; cumulNucFeats{i,3}]; %#ok<AGROW>
    groupsnucs = [groupsnucs; repmat(yeastTypes(i), size(cumulNucFeats{i,1}))]; %#ok<AGROW>
end

h_bp = cell([7 1]);
for i = 1:size(h_bp,1)
    h_bp{i} = figure('Visible', 'off');
end

figure(h_bp{1}); boxplot(boxIntSum, groups, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Total spot intensities'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;
figure(h_bp{2}); boxplot(boxIntRel, groups, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Relative spot intensities'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;
figure(h_bp{3}); boxplot(boxVol, groups, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Spot volumes'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;
figure(h_bp{4}); boxplot(boxVolRel, groups, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Relative spot volumes'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;

figure(h_bp{5}); boxplot(boxIntNuc, groupsnucs, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Nucleus Intensities'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;
figure(h_bp{6}); boxplot(boxBgNuc, groupsnucs, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Nucleus background intensities'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;
figure(h_bp{7}); boxplot(boxVolNuc, groupsnucs, 'Colors', cmap_yeast, 'Notch', 'on', 'Jitter', 0.3); ylabel('Nucleus volumes'); drawnow;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),h(j).Color,'FaceAlpha',.5);
end
drawnow;

% Write boxplots
for i = 1:size(h_bp,1)
    set(h_bp{i}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
end
saveas(h_bp{1}, [graphDir, filesep, 'boxplot_spot_intensities.pdf']); close(h_bp{1});
saveas(h_bp{2}, [graphDir, filesep, 'boxplot_spot_relative_intensities.pdf']); close(h_bp{2});
saveas(h_bp{3}, [graphDir, filesep, 'boxplot_spot_volumes.pdf']); close(h_bp{3});
saveas(h_bp{4}, [graphDir, filesep, 'boxplot_spot_relative_volumes.pdf']); close(h_bp{4});

saveas(h_bp{5}, [graphDir, filesep, 'boxplot_nucleus_intensities.pdf']); close(h_bp{5});
saveas(h_bp{6}, [graphDir, filesep, 'boxplot_nucleus_background.pdf']); close(h_bp{6});
saveas(h_bp{7}, [graphDir, filesep, 'boxplot_nucleus_volumes.pdf']); close(h_bp{7});

% Grouped histogram of number of nuclei per strain and number of foci
hyeast = figure('Visible', 'off');
byeast = bar(0:nbFociMax+1, nbFociHist');
legend(yeastTypes, 'interpreter', 'none'); xlabel('Nb foci per nucleus'); ylabel('Count');

for i = 1:nbYeastTypes
    byeast(i).FaceColor = cmap_yeast(i,:);
    byeast(i).EdgeColor = 'k';
end
set(hyeast, 'Units','centimeters', 'InvertHardcopy', 'off', ...
    'Position',[25 5 19 15], 'PaperSize',[19 15],...
    'PaperPositionMode','auto', 'Color', [1 1 1]);
saveas(hyeast, [graphDir, filesep, 'histograms.pdf']); close(hyeast);

% Write histogram table
fhistid = fopen([graphDir, filesep, 'histograms.txt'], 'w');
fprintf(fhistid, 'Yeast');
for inbfoci = 1:nbFociMax+1
    fprintf(fhistid, '\t%d foci', inbfoci-1);
end
fprintf(fhistid, '\n');
for iy = 1:nbYeastTypes
    fprintf(fhistid, '%s', yeastTypes{iy});
    for inbfoci = 1:nbFociMax+1
        fprintf(fhistid, '\t%d', nbFociHist(iy, inbfoci));
    end
    if iy < nbYeastTypes
        fprintf(fhistid, '\n');
    end
end
fclose(fhistid);

% Grouped boxplots of intensities by strain and number of foci
[~, ~, ic] = unique(groups);
dummynb = cat(1,cumulSpotsNb{:});
for inb = 1:nbFociMax
    idx = dummynb == inb;
    dummygps(idx) = ic(idx) + (inb-1)*nbYeastTypes; %#ok<AGROW>
end

positions = [];
for inb = 1:nbFociMax
    positions = [positions, repmat(inb, 1, nbYeastTypes) + (0:0.33/(nbYeastTypes-1):0.33)]; %#ok<AGROW>
end

xticksposi = [];
for inb= 1:nbFociMax
    xticksposi = [xticksposi, mean(positions((inb-1)*nbYeastTypes+(1:nbYeastTypes)))]; %#ok<AGROW>
end

positions(setdiff(1:(nbFociMax*nbYeastTypes), dummygps)) = [];
colors = repmat(cmap_yeast, nbFociMax, 1);
colors(setdiff(1:(nbFociMax*nbYeastTypes), dummygps), :) = [];

% Check dummygps to make sure legend is ok
[uniquegps, idum] = unique(dummygps);
[~, orddum] = sort(idum);
yeast_calls = mod(uniquegps(orddum), nbYeastTypes)+1;
first_calls = zeros(nbYeastTypes, 1);
for iy = 1:nbYeastTypes
    first_calls(iy) = find(yeast_calls == iy, 1, 'first');
end

h_bp = cell([4 1]);

h_bp{1} = figure('Visible', 'off'); boxplot(boxIntRel,dummygps, 'positions', positions);pause(1);
ax = h_bp{1}.Children;
set(ax,'xtick', xticksposi); set(ax,'xticklabel', 1:nbFociMax);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   [~, jid] = min(abs(mean(h(j).XData) - positions));
   patch(get(h(j),'XData'),get(h(j),'YData'),colors(jid,:),'FaceAlpha',.5);
end
drawnow;pause(1);
c = get(ax, 'Children'); legend(c(first_calls), yeastTypes, 'interpreter', 'none'); xlabel('Nb foci per nucleus'); ylabel('Relative spot intensities');

h_bp{2} = figure('Visible', 'off'); boxplot(boxIntSum,dummygps, 'positions', positions);pause(1);
ax = h_bp{2}.Children;
set(ax,'xtick', xticksposi); set(ax,'xticklabel', 1:nbFociMax);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   [~, jid] = min(abs(mean(h(j).XData) - positions));
   patch(get(h(j),'XData'),get(h(j),'YData'),colors(jid,:),'FaceAlpha',.5);
end
drawnow;pause(1);
c = get(ax, 'Children'); legend(c(first_calls), yeastTypes, 'interpreter', 'none'); xlabel('Nb foci per nucleus'); ylabel('Total spot intensities');

h_bp{3} = figure('Visible', 'off'); boxplot(boxVol,dummygps, 'positions', positions);pause(1);
ax = h_bp{3}.Children;
set(ax,'xtick', xticksposi); set(ax,'xticklabel', 1:nbFociMax);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   [~, jid] = min(abs(mean(h(j).XData) - positions));
   patch(get(h(j),'XData'),get(h(j),'YData'),colors(jid,:),'FaceAlpha',.5);
end
drawnow;pause(1);
c = get(ax, 'Children'); legend(c(first_calls), yeastTypes, 'interpreter', 'none'); xlabel('Nb foci per nucleus'); ylabel('Spot volumes');

h_bp{4} = figure('Visible', 'off'); boxplot(boxVolRel,dummygps, 'positions', positions);pause(1);
ax = h_bp{4}.Children;
set(ax,'xtick', xticksposi); set(ax,'xticklabel', 1:nbFociMax);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   [~, jid] = min(abs(mean(h(j).XData) - positions));
   patch(get(h(j),'XData'),get(h(j),'YData'),colors(jid,:),'FaceAlpha',.5);
end
drawnow;pause(1);
c = get(ax, 'Children'); legend(c(first_calls), yeastTypes, 'interpreter', 'none'); xlabel('Nb foci per nucleus'); ylabel('Relative spot volumes');

% Write grouped boxplots
for i = 1:size(h_bp,1)
    set(h_bp{i}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
end
saveas(h_bp{2}, [graphDir, filesep, 'grouped_boxplot_spot_intensities.pdf']); close(h_bp{2});
saveas(h_bp{1}, [graphDir, filesep, 'grouped_boxplot_spot_relative_intensities.pdf']); close(h_bp{1});
saveas(h_bp{3}, [graphDir, filesep, 'grouped_boxplot_spot_volumes.pdf']); close(h_bp{3});
saveas(h_bp{4}, [graphDir, filesep, 'grouped_boxplot_spot_relative_volumes.pdf']); close(h_bp{4});

% Scatterhists
h_sc = cell([5 1]);
h_sc{1} = figure('Visible', 'off'); scatterhist(boxIntRel(boxNbFoci == 1), boxD2CRel(boxNbFoci == 1), 'Group', groups(boxNbFoci == 1), 'Marker', '.', 'Kernel', 'off', 'Color', cmap_yeast);
xlabel('Relative intensity');
ylabel('Relative distance to center');

h_sc{2} = figure('Visible', 'off'); scatterhist(boxIntRel(boxNbFoci == 1), boxD2CRel_norad(boxNbFoci == 1), 'Group', groups(boxNbFoci == 1), 'Marker', '.', 'Kernel', 'off', 'Color', cmap_yeast);
xlabel('Relative intensity');
ylabel('Relative distance to center (no radius)');

h_sc{3} = figure('Visible', 'off'); scatterhist(boxIntRel(boxNbFoci == 1), boxD2C(boxNbFoci == 1), 'Group', groups(boxNbFoci == 1), 'Marker', '.', 'Kernel', 'off', 'Color', cmap_yeast);
xlabel('Relative intensity');
ylabel('Distance to center µm');

h_sc{4} = figure('Visible', 'off'); scatterhist(boxIntRel(boxNbFoci == 1), boxD2B(boxNbFoci == 1), 'Group', groups(boxNbFoci == 1), 'Marker', '.', 'Kernel', 'off', 'Color', cmap_yeast);
xlabel('Relative intensity');
ylabel('Distance to border µm');

h_sc{5} = figure('Visible', 'off'); scatterhist(boxIntRel(boxNbFoci == 1), boxD2B_norad(boxNbFoci == 1), 'Group', groups(boxNbFoci == 1), 'Marker', '.', 'Kernel', 'off', 'Color', cmap_yeast);
xlabel('Relative intensity');
ylabel('Distance to border (no radius) µm');

% Write scatterhists
for i = 1:size(h_sc,1)
    set(h_sc{i}, 'Units','centimeters', 'InvertHardcopy', 'off', ...
        'Position',[25 5 19 15], 'PaperSize',[19 15],...
        'PaperPositionMode','auto', 'Color', [1 1 1]);
end
saveas(h_sc{1}, [graphDir, filesep, 'scatterhist_relInt_relDistCenter.pdf']); close(h_sc{1});
saveas(h_sc{2}, [graphDir, filesep, 'scatterhist_relInt_relDistCenter_norad.pdf']); close(h_sc{2});
saveas(h_sc{3}, [graphDir, filesep, 'scatterhist_relInt_DistCenter.pdf']); close(h_sc{3});
saveas(h_sc{4}, [graphDir, filesep, 'scatterhist_relInt_DistBorder.pdf']); close(h_sc{4});
saveas(h_sc{5}, [graphDir, filesep, 'scatterhist_relInt_DistBorder_norad.pdf']); close(h_sc{5});

end
