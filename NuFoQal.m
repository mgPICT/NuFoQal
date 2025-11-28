function NuFoQal()

clear variables

filepath = fileparts(mfilename('fullpath'));
cd(filepath);

addpath(filepath);
addpath(genpath([filepath, filesep,'utilities']));

%% ================================================= %%
% User interface default colors
colorFgd=1-[0 0 0];
colorBgd=1-[0.85 0.85 0.87];
colorBut=1-[1 1 1];
colorCurie = [1 0.5216 0];


cprintf(sprintf('*[%.4f,%.4f,%.4f]',colorCurie(1), colorCurie(2), colorCurie(3)), '============= NuFoQal / Institut Curie - UMR 3664 =============\n');

%% ================================================= %%

% Create and hide the GUI as it is being constructed.
frontpanel = figure('Visible','on','Position',[00,000,500,400],...
    'Color',colorBgd,'Resize','off',...
    'Name', 'NuFoQal - Institut Curie',...  % Title figure
    'NumberTitle', 'off',... % Do not show figure number
    'MenuBar', 'none');
movegui(frontpanel, 'center');

START_PATH=filepath;
ANALYSE_FOLDER=0;
ANALYZE_PATHNAME='';
ANALYZE_FILENAME='';
RESULTS_PATHNAME='';
RESULTS_FILENAME='';
SAVE_IMGS=false;
TWO_SETS=false;
USE_NUCS=false;
PARAMS= {...
        '0.065';... resXY
        '0.200';... resZ
        '62';...    max nucleus diameter
        '2:3';...   spot size
        '0';...     Z slices to skip
        'GFP';...   signal tag
        'false';... low signal
        'true';...  use cellpose
        '';...      nb of std above mean
        '0';...     nb of pixels to add to contours
        'yAT';...   yeast indicator
        'DAPI';...  nuclei tag
        'cyto';...  model
        '30';...    cellDiameter
        'false'}; % useEnsemble

setappdata(gcf,'START_PATH',START_PATH);
setappdata(gcf,'ANALYSE_FOLDER',ANALYSE_FOLDER);
setappdata(gcf,'ANALYZE_PATHNAME',ANALYZE_PATHNAME);
setappdata(gcf,'ANALYZE_FILENAME',ANALYZE_FILENAME);
setappdata(gcf,'RESULTS_PATHNAME',RESULTS_PATHNAME);
setappdata(gcf,'RESULTS_FILENAME',RESULTS_FILENAME);
setappdata(gcf,'SAVE_IMGS', SAVE_IMGS);
setappdata(gcf,'TWO_SETS',TWO_SETS);
setappdata(gcf,'USE_NUCS', USE_NUCS);
setappdata(gcf,'PARAMS',PARAMS);
setappdata(gcf,'COLORS', cat(1, colorFgd, colorBgd, colorBut, colorCurie));

%% ================================================= %%
% The command line below defined the spatial structure of the user
% interface. Their related functions are defined at the end of this file.

%% ================================================= %%
%
%                    Panel NuFoQal
%
%% ================================================= %%

%% ================================================= %%

%% ================================================= %%
LoadData_panel = uipanel('Parent',frontpanel,'Title','Paths and data parameters',...
    'Position',[.05 .05 .5 .85], 'FontWeight', 'bold', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);
LoadData_File_or_Folder_text = uicontrol('Parent',frontpanel,'Style','text',...
    'String','Nuclear Foci Quantification and localization','Visible','on',...
    'FontSize',17', 'FontWeight', 'Bold',...
    'Units','normalized','Position',[0,0.9,1,0.1],...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorCurie); %#ok<NASGU> 

%% ================================================= %%
LoadData_File_or_Folder_buttonGroup = uibuttongroup('Parent',LoadData_panel,...
    'Visible','on',...
    'Units','normalized','Position',[0.05,0.45,0.9,0.5],...
    'SelectionChangeFcn',@LoadData_File_or_Folder_buttonGroup_Callback, ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd);

% LoadData_File_or_Folder_text = uicontrol('Parent',LoadData_File_or_Folder_buttonGroup,'Style','text',...
%     'String','Parameters','Visible','on',...
%     'Units','normalized','Position',[0.05,0.72,0.90,0.25],...
%     'FontWeight', 'bold', ...
%     'BackgroundColor',colorBgd,...
%     'ForegroundColor',colorFgd);

% Create two radio buttons in the button group.
u1 = uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','Radio','String',' Single file analysis',...
    'Units','normalized','Position',[0.05,0.77,0.90,0.2],...
    'HandleVisibility','on', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd); %#ok<NASGU> 
u2 = uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','Radio','String',' Folder analysis',...
    'Units','normalized','Position',[0.05,0.6,0.90,0.2],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd); %#ok<NASGU>

% change appearance of radioButtons
pathrbon=strcat([filepath,filesep,'utilities',filesep,'rb_on.png']);%'institut_curie.png']);
[imrb_on, ~, alpha_on] = imread(pathrbon);
imu1 = uipanel('Parent',frontpanel,...
    'Position',[0.1,0.735,0.033,0.041],...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'Borderwidth', 0);
ax = axes('parent',imu1,'units','normalized','position',[0 0 1 1]);
f1=image(imrb_on,'parent',ax);
axis off;
set(f1, 'AlphaData', alpha_on);

pathrboff=strcat([filepath,filesep,'utilities',filesep,'rb_off.png']);%'institut_curie.png']);
[imrb_off, ~, alpha_off] = imread(pathrboff);
imu2 = uipanel('Parent',frontpanel,...
    'Position',[0.1,0.67,0.033,0.0415],...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'BorderWidth', 0);
ax = axes('parent',imu2,'units','normalized','position',[0 0 1 1]);
f2=image(imrb_off,'parent',ax);
axis off;
set(f2, 'AlphaData', alpha_off);

% Create checkboxes for data settings
c0 = uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','CheckBox','String',' Nuclei on other channel',...
    'Units','normalized','Position',[0.05,0.37,0.90,0.15],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'Callback',@checkbox_nucleichannel_callback); %#ok<NASGU> 
c0 = uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','CheckBox','String',' Use deconvolution images',...
    'Units','normalized','Position',[0.05,0.2,0.90,0.15],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'Callback',@checkbox_useTwoSets_callback); %#ok<NASGU> 
c1 = uicontrol('Parent',LoadData_File_or_Folder_buttonGroup, ...
    'Style','CheckBox','String',' Save segmented images',...
    'Units','normalized','Position',[0.05,0.03,0.90,0.15],...
    'HandleVisibility','off', ...
    'BackgroundColor',colorBgd,...
    'ForegroundColor',colorFgd,...
    'Callback',@checkbox_saveimgs_callback); %#ok<NASGU> 


%% ================================================= %%
Folder_pushbutton=uicontrol('Parent',LoadData_panel,'Style','pushbutton',...
    'String','Choose data location',...
    'Units','normalized','Position',[0.05,0.25,0.9,0.15],...
    'TooltipString', 'Folder containing the data or file to analysis',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@Folder_pushbutton_Callback}); %#ok<NASGU> 

%% ================================================= %%
Result_pushbutton=uicontrol('Parent',LoadData_panel,'Style','pushbutton',...
    'String','Choose results location',...
    'Units','normalized','Position',[0.05,0.05,0.9,0.15],...
    'TooltipString', 'Folder in which the analysis results will be written',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@Result_pushbutton_Callback}); %#ok<NASGU> 

%% ================================================= %%
% To easily find the location of the icon image without considering O/S and
% display it in the GUI
pathImgUI=strcat([filepath,filesep,'Logo_PICT_neg.png']);%'institut_curie.png']);
[imageLogo, ~, alpha] = imread(pathImgUI);
Logo_panel=uipanel('Parent',frontpanel,...
    'Position',[0.60,0.66,0.15,0.1388],...
    'BackgroundColor',colorBgd,...
    'BorderWidth', 0, ...
    'ForegroundColor',colorFgd);
ax = axes('parent',Logo_panel,'units','normalized','position',[0 0 1 1]);
fp=image(imageLogo,'parent',ax);
axis off;
set(fp, 'AlphaData', alpha,'Interpolation', 'bilinear');

pathImgUI=strcat([filepath,filesep,'LogosCurie_neg.png']);%'institut_curie.png']);
[imageLogo, ~, alpha] = imread(pathImgUI);
Logo_panel=uipanel('Parent',frontpanel,...
    'Position',[0.78,0.66,0.17,0.216335],...
    'BackgroundColor',colorBgd,...
    'BorderWidth', 0, ...
    'ForegroundColor',colorFgd);
ax = axes('parent',Logo_panel,'units','normalized','position',[0 0 1 1]);
fc=image(imageLogo,'parent',ax);
axis off;
set(fc, 'AlphaData', alpha,'Interpolation', 'bilinear');

%% ================================================= %%
Params_pushbutton = uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','Set parameters',...
    'Units','normalized','Position',[0.60,0.51,0.35,0.1],...
    'TooltipString', 'Open a window to enter the analysis parameters',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorFgd,...
    'Callback',{@Params_pushbutton_Callback}); %#ok<NASGU> 

%% ================================================= %%
RunAnalysis_pushbutton=uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','ANALYSIS', 'FontSize', 20, 'fontWeight', 'Bold', ...
    'Units','normalized','Position',[0.60,0.35,0.350,0.1215],...
    'TooltipString', 'Run the analysis on the selected files',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorCurie,...
    'Callback',{@RunAnalysis_pushbutton_Callback}); %#ok<NASGU> 

FilterResults_pushbutton=uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','CHECK', 'FontSize', 20, 'fontWeight', 'Bold', ...
    'Units','normalized','Position',[0.60,0.225,0.350,0.1215],...
    'TooltipString', 'Supervised result filtering on the selected files',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorCurie,...
    'Callback',{@FilterResults_pushbutton_Callback}); %#ok<NASGU> 

DrawPlots_pushbutton=uicontrol('Parent',frontpanel,'Style','pushbutton',...
    'String','GRAPHS', 'FontSize', 20, 'fontWeight', 'Bold', ...
    'Units','normalized','Position',[0.60,0.095,0.350,0.1215],...
    'TooltipString', 'Generate figures from the analysis conducted on the selected files',...
    'BackgroundColor',colorBut,...
    'ForegroundColor',colorCurie,...
    'Callback',{@DrawPlots_pushbutton_Callback}); %#ok<NASGU> 
% End of the qFOCIpanel
%% ===========  End of the qFOCIpanel   ============ %%
%%================================================== %%
%%================================================== %%
%%================================================== %%

end%function

%% function LoadData_File_or_Folder_buttonGroup_Callback(source,eventdata)
% Set the king of analysis: only on a single file (ANALYSE_FOLDER=0) or
% on a whole folder (ANALYSE_FOLDER=1)
function LoadData_File_or_Folder_buttonGroup_Callback(source, ~)
    str=get(get(source,'SelectedObject'),'String');
    
    ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
    if strcmp(str,' Single file analysis')
        ANALYSE_FOLDER=0;
        disp(' Single file analysis')
    end
    
    if strcmp(str,' Folder analysis')
        ANALYSE_FOLDER=1;
        disp(' Folder analysis')
    end
    
    % Invert radiobutton images
    rb_single = gcbf().Children(7);
    rb_folder = gcbf().Children(8);
    
    tmp_im = rb_single.Children.Children.CData;
    tmp_al = rb_single.Children.Children.AlphaData;
    rb_single.Children.Children.CData = rb_folder.Children.Children.CData;
    rb_single.Children.Children.AlphaData = rb_folder.Children.Children.AlphaData;
    rb_folder.Children.Children.CData = tmp_im;
    rb_folder.Children.Children.AlphaData = tmp_al;
    
    clear tmp_im tmp_al
    
    setappdata(gcbf,'ANALYSE_FOLDER',ANALYSE_FOLDER);
end

%% function checkbox_nucleichannel_callback(source,eventdata)
% Set wether the segmentation is done using two sets of images or not
function checkbox_nucleichannel_callback(source, ~)
    val = get(source,'Value');
    
    if val
        cprintf('comment', 'Nuclei segmentation will be done on a specific channel: %s\n');
    else
        cprintf('comment', 'Nuclei segmentation will rely on nuclear background of spots channel.\n');
    end
    
    setappdata(gcbf, 'USE_NUCS', val);
end


%% function checkbox_useTwoSets_callback(source,eventdata)
% Set wether the segmentation is done using two sets of images or not
function checkbox_useTwoSets_callback(source, ~)
    val = get(source,'Value');

    TWO_SETS = val;
    if val
        cprintf('comment', 'Spots segmentation will be done with two different sets of images.\n');
    else
        cprintf('comment', 'Spots segmentations and quantifications will be done on a single set of images.\n');
    end
    
    setappdata(gcbf, 'TWO_SETS', TWO_SETS);
end


%% function checkbox_saveimgs_callback(source,eventdata)
% Set the saving of images: checked (true) unchecked (false)
function checkbox_saveimgs_callback(source, ~)
    val=get(source,'Value');
    
    SAVE_IMGS = val;
    if val
        cprintf('comment', ' Segmented images will be written in the result location.\n')
    else
        cprintf('comment', ' Segmented images will not be written.\n')
    end
    
    setappdata(gcbf,'SAVE_IMGS',SAVE_IMGS);
end

%% function Folder_pushbutton_Callback(hObject, eventdata, handles)
% Open a dialog box to specify where the files to analyzed are located. 
% Depending the kind of analysis (single file or all file in a folder), the
% program ask to specify the image or the folder to be analyzed.
function Folder_pushbutton_Callback(~, ~, ~)
    ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
    START_PATH=getappdata(gcbf,'START_PATH');
    ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
    ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
    TWO_SETS = getappdata(gcbf, 'TWO_SETS');
    ext = {'*.tif','TIF image (.tif)';'*.TIF','TIF image (.TIF)'};
    
    if isempty(ANALYZE_PATHNAME)
        cd(START_PATH);
    else
        if iscell(ANALYZE_PATHNAME)
            if isfolder(ANALYZE_PATHNAME{1})
                cd(ANALYZE_PATHNAME{1});
            else
                cd(START_PATH);
            end
        else
            if isfolder(ANALYZE_PATHNAME)
                cd(ANALYZE_PATHNAME);
            else
                cd(START_PATH);
            end
        end
    end
    
    if TWO_SETS
        if (ANALYSE_FOLDER)
            [ANALYZE_PATHNAME] = uigetdir('./','Select the folder containing the quantification set of images.');
            cd(ANALYZE_PATHNAME);
            [ANALYZE_PATHNAME2] = uigetdir(ANALYZE_PATHNAME,'Select the folder containing the spots segmentation set of images.');
        else
            [ANALYZE_FILENAME,ANALYZE_PATHNAME] = uigetfile(ext,'Select the quantification image file to analyse.');
            cd(ANALYZE_PATHNAME);
            [ANALYZE_FILENAME2,ANALYZE_PATHNAME2] = uigetfile(ext,'Select the spots segmentation image file to analyse.');
            ANALYZE_FILENAME = {ANALYZE_FILENAME, ANALYZE_FILENAME2};
        end
        ANALYZE_PATHNAME = {ANALYZE_PATHNAME, ANALYZE_PATHNAME2};
    else
        if (ANALYSE_FOLDER)
            [ANALYZE_PATHNAME] = uigetdir('./','Select the folder containing the images.');
        else
            [ANALYZE_FILENAME,ANALYZE_PATHNAME] = uigetfile(ext,'Select the image file to analyse.');
        end
    end
    
    cd(START_PATH);
    setappdata(gcbf,'ANALYZE_PATHNAME',ANALYZE_PATHNAME);
    setappdata(gcbf,'ANALYZE_FILENAME',ANALYZE_FILENAME);
end

function Result_pushbutton_Callback(~, ~, ~)
    START_PATH=getappdata(gcbf,'START_PATH');
    ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
    RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
    
    if isempty(RESULTS_PATHNAME)
        if isempty(ANALYZE_PATHNAME)
            cd(START_PATH);
        else
            if iscell(ANALYZE_PATHNAME)
                if isfolder(ANALYZE_PATHNAME{1})
                    cd(ANALYZE_PATHNAME{1});
                else
                    cd(START_PATH);
                end
            else
                if isfolder(ANALYZE_PATHNAME)
                    cd(ANALYZE_PATHNAME);
                else
                    cd(START_PATH);
                end
            end
        end
    else
        if isfolder(RESULTS_PATHNAME)
            cd(RESULTS_PATHNAME);
        else
            cd(START_PATH);
        end
    end
    
    [RESULTS_PATHNAME] = uigetdir('./','Select (or create) the output folder.');
    
    cd(START_PATH);
    setappdata(gcbf,'RESULTS_PATHNAME',RESULTS_PATHNAME);
end

%% function Params_pushbutton_Callback(hObject, eventdata, handles)
% Open a dialog box to set the different wavelengths to re-align
function Params_pushbutton_Callback(~, ~, ~)
    PARAMS=getappdata(gcbf,'PARAMS');
    COLORS = getappdata(gcbf, 'COLORS');
    gui_params(gcbf, PARAMS,COLORS(1,:), COLORS(2,:), COLORS(3,:), COLORS(4,:));
end


%% function RunAnalysis_pushbutton_Callback(hObject, eventdata, handles)
% Run the qFOCI analysis
function RunAnalysis_pushbutton_Callback(~, ~, ~)
    START_PATH = getappdata(gcbf,'START_PATH');
    ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
    ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
    ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
    RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
    SAVE_IMGS = getappdata(gcbf,'SAVE_IMGS');
    PARAMS=getappdata(gcbf,'PARAMS');
    USE_NUCS=getappdata(gcbf,'USE_NUCS');
    
    cd(START_PATH);
    error = nufoqal_analysis(ANALYZE_PATHNAME, RESULTS_PATHNAME, SAVE_IMGS, PARAMS, ~ANALYSE_FOLDER, ANALYZE_FILENAME, USE_NUCS);
    if error > 0
        fprintf("NuFoQal analysis returned with error: #%d", error);
    end
end

function FilterResults_pushbutton_Callback(~, ~, ~)
    START_PATH = getappdata(gcbf,'START_PATH');
    ANALYSE_FOLDER=getappdata(gcbf,'ANALYSE_FOLDER');
    ANALYZE_PATHNAME=getappdata(gcbf,'ANALYZE_PATHNAME');
    ANALYZE_FILENAME=getappdata(gcbf,'ANALYZE_FILENAME');
    RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
    SAVE_IMGS = getappdata(gcbf,'SAVE_IMGS');
    USE_NUCS = getappdata(gcbf, 'USE_NUCS');
    PARAMS=getappdata(gcbf,'PARAMS');
    COLORS = getappdata(gcbf, 'COLORS');
    
    cd(START_PATH);
    error = supervised_filtering_results_mod(ANALYZE_PATHNAME, RESULTS_PATHNAME, SAVE_IMGS, PARAMS, ~ANALYSE_FOLDER, ANALYZE_FILENAME, START_PATH, USE_NUCS, COLORS(1,:), COLORS(2,:), COLORS(3,:), COLORS(4,:));
    if error > 0
        fprintf("Supervised filtering results returned with error: %s", error);
    end
end

function DrawPlots_pushbutton_Callback(~, ~, ~)
    START_PATH = getappdata(gcbf,'START_PATH');
    RESULTS_PATHNAME=getappdata(gcbf,'RESULTS_PATHNAME');
    PARAMS=getappdata(gcbf,'PARAMS');
    
    cd(START_PATH);
    error = generate_figures(RESULTS_PATHNAME, PARAMS);
    if error > 0
        fprintf("Generate figures returned with error: %s", error);
    end
end
