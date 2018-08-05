%% Function: fw_gui
%%
%% Description: GUI to access image-space functionality of fat-water toolbox
%%
%% Author: Diego Hernando
%% Date created: November 20, 2011
%% Date last modified: February 10, 2012

function fw_gui

global imDataParams;
global outParams;
global butProcess;
global butSave;
global refAlgoPopup;
global butProcessRef;
global butSaveRef;
global butPostProcess;
global butSavePost;
global postAlgoPopup;
global imwait waitAxes


initGlobal();
%   Set Number of tabs and tab labels.  Make sure the number of tab labels
%   match the HumberOfTabs setting.
NumberOfTabs = 4;               % Number of tabs to be generated
TabLabels = {'Load data'; 'Processing'; 'Refinements'; 'Postprocessing'};
if size(TabLabels,1) ~= NumberOfTabs
  errordlg('Number of tabs and tab labels must be the same','Setup Error');
  return
end

%   Get user screen size
set(0,'Units','pixels');
SC = get(0, 'ScreenSize');
MaxMonitorX = SC(3);
MaxMonitorY = SC(4);

%   Set the figure window size values
MainFigScale = .8;          % Change this value to adjust the figure size
MaxWindowX = min(round(MaxMonitorX*MainFigScale),1200);
MaxWindowY = min(round(MaxMonitorY*MainFigScale),750);
XBorder = min(40,(MaxMonitorX-MaxWindowX)/2);
YBorder = (MaxMonitorY-MaxWindowY)/2; 
TabOffset = 0;              % This value offsets the tabs inside the figure.
ButtonHeight = 40;
PanelWidth = MaxWindowX-2*TabOffset+4;
PanelHeight = MaxWindowY-ButtonHeight-2*TabOffset;
ButtonWidth = round((PanelWidth-NumberOfTabs)/NumberOfTabs);

%   Set the color varables.  
White = [1  1  1];            % White - Selected tab color     
BGColor = .8*White;           % Light Grey - Background color

%%   Create a figure for the tabs
hTabFig = figure(...
    'Units', 'pixels',...
    'Toolbar', 'none',...
    'Position',[ XBorder, YBorder, MaxWindowX, MaxWindowY ],...
    'NumberTitle', 'off',...
    'Name', 'Fat-Water Toolbox',...
    'MenuBar', 'none',...
    'Resize', 'off',...
    'DockControls', 'off',...
    'Color', White);

%set(hTabFig,'Toolbar','figure');

%%   Define a cell array for panel and pushbutton handles, pushbuttons labels and other data
%   rows are for each tab + two additional rows for other data
%   columns are uipanel handles, selection pushbutton handles, and tab label strings - 3 columns.
TabHandles = cell(NumberOfTabs,3);
TabHandles(:,3) = TabLabels(:,1);
%   Add additional rows for other data
TabHandles{NumberOfTabs+1,1} = hTabFig;         % Main figure handle
TabHandles{NumberOfTabs+1,2} = PanelWidth;      % Width of tab panel
TabHandles{NumberOfTabs+1,3} = PanelHeight;     % Height of tab panel
TabHandles{NumberOfTabs+2,1} = 0;               % Handle to default tab 2 content(set later)
TabHandles{NumberOfTabs+2,2} = White;           % Selected tab Color
TabHandles{NumberOfTabs+2,3} = BGColor;         % Background color

%%   Build the Tabs
for TabNumber = 1:NumberOfTabs
  
  % create a UIPanel   
  TabHandles{TabNumber,1} = uipanel('Units', 'pixels', ...
                                    'Visible', 'off', ...
                                    'Backgroundcolor', White, ...
                                    'BorderWidth',1, ...
                                    'Position', [TabOffset TabOffset ...
                      PanelWidth PanelHeight]);

  % create a selection pushbutton
  TabHandles{TabNumber,2} = uicontrol('Style', 'pushbutton',...
                                      'Units', 'pixels', ...
                                      'BackgroundColor', BGColor, ...
                                      'Position', [TabOffset+(TabNumber-1)*ButtonWidth PanelHeight+TabOffset...
                      ButtonWidth ButtonHeight], ...          
                                      'String', TabHandles{TabNumber,3},...
                                      'HorizontalAlignment', 'center',...
                                      'FontName', 'arial',...
                                      'FontWeight', 'bold',...
                                      'FontSize', 10);

end

%%   Define the callbacks for the Tab Buttons
%   All callbacks go to the same function with the additional argument being the Tab number
for CountTabs = 1:NumberOfTabs
  set(TabHandles{CountTabs,2}, 'callback', ...
                    {@TabSelectCallback, CountTabs});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
uicontrol('Parent', TabHandles{1,1}, ...
          'Units', 'pixels', ...
          'Position', [round(PanelWidth/2)-100 PanelHeight-2*ButtonHeight-50 200 ButtonHeight], ...
          'String', 'Load data', ...
          'Callback', @OpenImageCallback , ...
          'Style', 'pushbutton',...
          'HorizontalAlignment', 'center',...
          'FontName', 'arial',...
          'FontWeight', 'bold',...
          'FontSize', 12);


imwait = imread('hour1.png');
waitAxes(1) = axes('Parent', TabHandles{1,1}, ...
                 'Units', 'pixels', ...
                 'Visible','off',...
                 'Position', [PanelWidth-100 PanelHeight-100 40 60]);

waitAxes(2) = axes('Parent', TabHandles{2,1}, ...
                 'Units', 'pixels', ...
                 'Visible','off',...
                 'Position', [PanelWidth-100 PanelHeight-100 40 60]);

waitAxes(3) = axes('Parent', TabHandles{3,1}, ...
                 'Units', 'pixels', ...
                 'Visible','off',...
                 'Position', [PanelWidth-100 PanelHeight-100 40 60]);

waitAxes(4) = axes('Parent', TabHandles{4,1}, ...
                 'Units', 'pixels', ...
                 'Visible','off',...
                 'Position', [PanelWidth-100 PanelHeight-100 40 60]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Display it - Put the handle in TabHandles so that it can be deleted later 
uicontrol('Style', 'popup',...
          'Parent', TabHandles{2,1}, ...
          'String', '2 point, flexible TEs (Berglund)|3 point, multi-seed (Berglund)|3+ point, multipeak (Tsao/Jiang)|3 point, multi-resol (Tsao/Jiang)|3+ point, graphcut (Hernando)|3 point, golden section (Lu)|3+ point, multipeak RG (Hu)|3+ point, multipeak (Sharma)',...
          'Position', [20 PanelHeight-120 280 50],...
          'FontSize',12,...
          'Callback', @setAlgo); 



%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butProcess = uicontrol('Parent', TabHandles{2,1}, ...
                       'Units', 'pixels', ...
                       'Position', [20 PanelHeight-60 280 40], ...
                       'String', 'Process Data', ...
                       'Callback', @processData , ...
                       'Style', 'pushbutton',...
                       'HorizontalAlignment', 'center',...
                       'FontName', 'arial',...
                       'Enable', 'off',...
                       'FontWeight', 'bold',...
                       'FontSize', 12);


%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butSave = uicontrol('Parent', TabHandles{2,1}, ...
                    'Units', 'pixels', ...
                    'Position', [PanelWidth-260 20 200 40], ...
                    'String', 'Save Results', ...
                    'Callback', @saveResults , ...
                    'Style', 'pushbutton',...
                    'HorizontalAlignment', 'center',...
                    'FontName', 'arial',...
                    'Enable', 'off',...
                    'FontWeight', 'bold',...
                    'FontSize', 12);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Define Tab 3 content


%   Display it - Put the handle in TabHandles so that it can be deleted later 
refAlgoPopup = uicontrol('Style', 'popup',...
                         'Parent', TabHandles{3,1}, ...
                         'String', 'Mixed magnitude/complex fitting |Noncartesian deblurring',...
                         'Position', [20 PanelHeight-120 280 50],...
                         'FontSize',12,...
                         'Callback', @setAlgoRef); 



%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butProcessRef = uicontrol('Parent', TabHandles{3,1}, ...
                          'Units', 'pixels', ...
                          'Position', [20 PanelHeight-60 280 40], ...
                          'String', 'Process Data', ...
                          'Style', 'pushbutton',...
                          'HorizontalAlignment', 'center',...
                          'FontName', 'arial',...
                          'Enable', 'off',...
                          'FontWeight', 'bold',...
                          'Callback', @processDataRef , ...
                          'FontSize', 12);
 

%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butSaveRef = uicontrol('Parent', TabHandles{3,1}, ...
                    'Units', 'pixels', ...
                    'Position', [PanelWidth-260 20 200 40], ...
                    'String', 'Save Results', ...
                    'Style', 'pushbutton',...
                    'HorizontalAlignment', 'center',...
                    'FontName', 'arial',...
                    'Enable', 'off',...
                    'FontWeight', 'bold',...
                    'Callback', @saveResultsRef , ...
                    'FontSize', 12);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Define Tab 4 content

%   Display it - Put the handle in TabHandles so that it can be deleted later 
postAlgoPopup = uicontrol('Style', 'popup',...
                          'Parent', TabHandles{4,1}, ...
                          'String', 'Compute NSA map|Compute fat fraction map',...
                          'Position', [20 PanelHeight-120 280 50],...
                          'Callback', @setAlgoPost,...
                          'FontSize',12); 



%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butPostProcess = uicontrol('Parent', TabHandles{4,1}, ...
                           'Units', 'pixels', ...
                           'Position', [20 PanelHeight-60 280 40], ...
                           'String', 'Post-Process Data', ...
                           'Callback', @postProcess , ...
                           'Style', 'pushbutton',...
                           'HorizontalAlignment', 'center',...
                           'FontName', 'arial',...
                           'Enable', 'off',...
                           'FontWeight', 'bold',...
                           'FontSize', 12);


%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butSavePost = uicontrol('Parent', TabHandles{4,1}, ...
                        'Units', 'pixels', ...
                        'Position', [PanelWidth-260 20 200 40], ...
                        'String', 'Save Results', ...
                        'Callback', @savePostResults , ...
                        'Style', 'pushbutton',...
                        'HorizontalAlignment', 'center',...
                        'FontName', 'arial',...
                        'Enable', 'off',...
                        'FontWeight', 'bold',...
                        'FontSize', 12);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Done  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Save the TabHandles in guidata
guidata(hTabFig,TabHandles);

%%   Make Tab 1 active
TabSelectCallback(0,0,1);

%%%%% WARNING %%%%%%
msgbox('This software is for research/development purposes. Not for clinical use.','Research software','warn')

end



function startWait( tabnum )
global imwait waitAxes
image(imwait,'Parent', waitAxes(tabnum));axis equal tight off;% freezeColors;
set(waitAxes(tabnum),'Visible','off')
drawnow;
end

function stopWait( tabnum )
global imwait waitAxes
image(255 + 0*imwait,'Parent', waitAxes(tabnum));axis equal tight off;% freezeColors;
set(waitAxes(tabnum),'Visible','off')
drawnow;
end


%%   Callback for Tab Selection
function TabSelectCallback(ev,source,SelectedTab)
%   All tab selection pushbuttons are greyed out and uipanels are set to
%   visible off, then the selected panel is made visible and it's selection
%   pushbutton is highlighted.

%   Set up some varables
TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
White = TabHandles{NumberOfTabs+2,2};            % White      
BGColor = TabHandles{NumberOfTabs+2,3};          % Light Grey

%   Turn all tabs off
for TabCount = 1:NumberOfTabs
  set(TabHandles{TabCount,1}, 'Visible', 'off');
  set(TabHandles{TabCount,2}, 'BackgroundColor', BGColor);
end

%   Enable the selected tab
set(TabHandles{SelectedTab,1}, 'Visible', 'on');        
set(TabHandles{SelectedTab,2}, 'BackgroundColor', White);

end


%%   Save water/fat separation results
function saveResults(ev,source)

global imDataParams algoParams outParams algoName wfAxes PicNameWithTag

uisave('outParams',['outParams_'  PicNameWithTag]);
end

%%   Save refined results
function saveResultsRef(ev,source)

global imDataParams algoParams outParamsRef algoName wfAxes PicNameWithTag

uisave('outParamsRef',['outParamsRef_'  PicNameWithTag]);
end


%% Save postprocessing results
function savePostResults(ev,source)

global imDataParams algoParams outParams algoName wfAxes PicNameWithTag postAlgoPopup butSavePost ff outNSA


ppChoice = get(postAlgoPopup,'Value');

switch ppChoice
  
 case 1
  
  uisave('outNSA',['outNSA_'  PicNameWithTag]);
  
 case 2

  uisave('ff',['fatFraction_'  PicNameWithTag]);
  
end
end

%%   process Dataset
function processData(ev,source)

global imDataParams algoParams outParams algoName wfAxes butSave butPostProcess butProcessRef
global outParamsRef wfAxesRef
global butSave butProcessRef butSaveRef butSavePost butPostProcess
global wfAxesRef postAxes outParamsRef


%   Get TabHandles from guidata and set some varables
TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};

%%%% START WAIT %%%%%
startWait( 2 );
%%%%%%%%%%%%%%%%%%%%%

eval(['outParams = ' algoName '(imDataParams,algoParams);']);

%%%% STOP WAIT %%%%%
stopWait( 2 );
%%%%%%%%%%%%%%%%%%%%%

if ~isempty(outParams)
  
  %   Set the axes and display the image    
  ImgSize = round(0.45*PanelWidth);
  %ImgSize = 700;
  dx = round(PanelWidth-ImgSize - 20);
  dy = round(PanelHeight/2-ImgSize/2);;
  
  try 
    delete(wfAxes);
  end
  wfAxes = axes('Parent', TabHandles{2,1}, ...
                'Units', 'pixels', ...
                'Position', [dx dy ImgSize ImgSize]);
  
  try 
    imwf = abs(cat(3,outParams.species(1).amps(:,:,ceil(end/2)),outParams.species(2).amps(:,:,ceil(end/2))));
  catch
    imwf = abs(cat(3,outParams.water(:,:,ceil(end/2)), outParams.fat(:,:,ceil(end/2))));  
  end
  
  sc = max(imwf(:));
  
  imagesc([imwf(:,:,1) sc*ones(size(imwf,1),4) imwf(:,:,2)],'Parent', wfAxes,[0 0.6*sc]);
  axis equal tight off;colormap gray; freezeColors;  
  % imagescn2(imwf,[0 0.75*sc],[],wfAxes,[dx dy ImgSize ImgSize]);axis equal tight off;colormap gray  

  fs = 24;
  title(wfAxes,'Water/Fat images','FontSize',fs)
  

  %% Clear refinements/postprocessing figures and data, update save buttons
  clear global outParamsRef
  set(butSave,'enable','on');
  set(butProcessRef,'enable','on');
  set(butSaveRef,'enable','off');
  set(butSavePost,'enable','off');
  set(butPostProcess,'enable','on');
  try 
    delete(wfAxesRef)
  end
  try 
    delete(postAxes)
  end

  
  

end

end

%%   process Dataset -- Refinements
function processDataRef(ev,source)

global imDataParams algoParamsRef outParams outParamsRef algoNameRef wfAxesRef butSaveRef butPostProcess refAlgoPopup



%   Get TabHandles from guidata and set some varables
TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};


%%%% START WAIT %%%%%
startWait( 3 );
%%%%%%%%%%%%%%%%%%%%%

try
  algoParamsRef.fieldmap = double(outParams.fieldmap);
catch 
  algoParamsRef.fieldmap = 10 + zeros(size(outParams.species(1).amps));
  disp('No fieldmap found, assuming zero fieldmap');
end

try
  algoParamsRef.r2starmap = double(outParams.r2starmap);
catch 
  algoParamsRef.r2starmap = zeros(size(outParams.species(1).amps));
  disp('No R2* map found, assuming all zeros');
end

ppChoice = get(refAlgoPopup,'Value');

if ppChoice==2 % Eggers' spiral deblurring
  

  if isfield(outParams,'fieldmap')
    
    is                    =  size( outParams.species(1).amps, 1 );
    
    lambda                =  algoParamsRef.lambda;        % Shape parameter
    delay                 =  algoParamsRef.delay;
    
    trajectory =  compute_constant_density_spiral_trajectory( algoParamsRef.nv, algoParamsRef.ns, is, lambda, delay );
    
    % Compute sampling density compensation weights
    density_compensation =  compute_constant_density_spiral_density( algoParamsRef.nv, algoParamsRef.ns, 2 * ceil( algoParamsRef.alpha * is / 2 ), is, lambda, delay, trajectory );
    
    trajectory =  is * trajectory + is / 2;
    
    algoParamsRef.trajectory           =  trajectory; 
    algoParamsRef.density_compensation =  density_compensation;
    
    gamma            = 42.58;
    [maxAmp,iMaxAmp] = max(algoParamsRef.species(2).relAmps);
    chemical_shift   = algoParamsRef.species(2).frequency(iMaxAmp); % Fat modeled as single peak here to keep it simple
    frequency_offset = gamma * chemical_shift * imDataParams.FieldStrength;
    
    outParams.species(1).fieldmap = - algoParamsRef.fieldmap;
    outParams.species(2).fieldmap = - algoParamsRef.fieldmap - frequency_offset;
    
    
    eval(['outParamsRef = ' algoNameRef '(outParams,algoParamsRef);']);
    
  else
    outParamsRef = [];
    disp('No fieldmap provided by processing algorithm (tab 2). Please use a different algorithm.');
  end
  
else % Hernando's mixed fitting. 
  
  eval(['outParamsRef = ' algoNameRef '(imDataParams,algoParamsRef);']);


end


%%%% STOP WAIT %%%%%
stopWait( 3 );
%%%%%%%%%%%%%%%%%%%%%

if ~isempty(outParamsRef)
  
  %   Set the axes and display the image    
  ImgSize = round(0.45*PanelWidth);
  %ImgSize = 700;
  dx = round(PanelWidth-ImgSize - 20);
  dy = round(PanelHeight/2-ImgSize/2);;
  
  try 
    delete(wfAxesRef);
  end
  wfAxesRef = axes('Parent', TabHandles{3,1}, ...
                'Units', 'pixels', ...
                'Position', [dx dy ImgSize ImgSize]);
  
  try 
    imwf = abs(cat(3,outParamsRef.species(1).amps(:,:,ceil(end/2)),outParamsRef.species(2).amps(:,:,ceil(end/2))));
  catch
    imwf = abs(cat(3,outParamsRef.water(:,:,ceil(end/2)), outParamsRef.fat(:,:,ceil(end/2))));  
  end
  
  sc = max(imwf(:));
  
  imagesc([imwf(:,:,1) sc*ones(size(imwf,1),4) imwf(:,:,2)],'Parent', wfAxesRef,[0 0.6*sc]);
  axis equal tight off;colormap gray; freezeColors;  
  % imagescn2(imwf,[0 0.75*sc],[],wfAxes,[dx dy ImgSize ImgSize]);axis equal tight off;colormap gray  

  fs = 24;
  title(wfAxesRef,'Water/Fat images','FontSize',fs)
  
  set(butSaveRef,'enable','on');
  set(butPostProcess,'enable','on');
end

end


%%   process Dataset
function postProcess(ev,source)

global imDataParams algoParams outParams outParamsRef algoName postAxes butSave postAlgoPopup butSavePost outNSA ff postParams allPostParams


%   Get TabHandles from guidata and set some varables
TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};


%%%% START WAIT %%%%%
startWait( 4 );
%%%%%%%%%%%%%%%%%%%%%

ppChoice = get(postAlgoPopup,'Value');

if ~exist('postParams')
  postParams = allPostParams(ppChoice).postParams
  
end


if exist('outParamsRef')
  if ~isempty(outParamsRef)
    curOutParams = outParamsRef;
  else
    curOutParams = outParams;  
  end
else
  curOutParams = outParams;  
end


switch ppChoice
  
 case 1
  
  curOutParams.include_r2star = postParams.include_r2star;
  
  outNSA = fw_nsa_map(imDataParams, algoParams, curOutParams);
  %ImgSize = 700;
  ImgSize = round(0.45*PanelWidth);
  dx = round(PanelWidth-ImgSize - 20);
  dy = round(PanelHeight/2-ImgSize/2);;
  
  try 
    delete(postAxes);
  end
  postAxes = axes('Parent', TabHandles{4,1}, ...
                  'Units', 'pixels', ...
                  'Position', [dx dy ImgSize ImgSize]);
  
  
  nsamap = cat(2,outNSA.species(1).nsa, outNSA.species(2).nsa);
  
  imagesc(nsamap(:,:,ceil(end/2)),[0 max(nsamap(:))]);axis equal tight off;colormap jet; freezeColors;  
  colorbar
  fs = 20;
  title(postAxes,'Number of Signal Averages (NSA) maps (Water/fat)','FontSize',fs);

  
  
  
 case 2

  ff = computeFF( curOutParams, postParams );

  ImgSize = round(0.6*PanelHeight);
  %  ImgSize = 500;
  dx = round(3*PanelWidth/4-ImgSize/2);
  dy = round(PanelHeight/2-ImgSize/2);;
  
  try 
    delete(postAxes);
  end
  postAxes = axes('Parent', TabHandles{4,1}, ...
                  'Units', 'pixels', ...
                  'Position', [dx dy ImgSize ImgSize]);
  
  imagesc(ff(:,:,ceil(end/2)),[0 100]);axis equal tight off;colormap gray;  freezeColors;  
  colorbar
  fs = 24;
  title(postAxes,'Fat fraction (%)','FontSize',fs);

end



%%%% STOP WAIT %%%%%
stopWait( 4 );
%%%%%%%%%%%%%%%%%%%%%



% $$$   
% $$$   eval(['outParams = ' algoName '(imDataParams,algoParams);']);
% $$$   
% $$$   %   Set the axes and display the image    


set(butSavePost,'enable','on');


end


%%   Open Image File Callback
function OpenImageCallback(ev,source)


global imDataParams butProcess PicNameWithTag
global butProcessRef butPostProcess wfAxes wfAxesRef postAxes outParams outParamsRef
global butSave butSaveRef butSavePost

%   Get TabHandles from guidata and set some varables
TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};


%   Two persistent varables are needed
persistent StartPicDirectory hImageAxes

%   Initilize the StartPicDirectory if first time through
if isempty(StartPicDirectory)
  StartPicDirectory = cd;
end

%   Get the file name from the user
[PicNameWithTag, PicDirectory] = uigetfile({'*.mat;*.MAT','Matlab Files'},...
                                           'Select a Matlab data file',StartPicDirectory);

if PicNameWithTag == 0,
  %   If User canceles then display error message
  errordlg('You should select a Matlab .mat File');
  return
end

%   Set the default directory to the currently selected directory
StartPicDirectory = PicDirectory;

%   Build path to file
PicFilePath = strcat(PicDirectory,PicNameWithTag);


%%%% START WAIT %%%%%
startWait( 1 );
%%%%%%%%%%%%%%%%%%%%%

%   Load the image
data1 = load(PicFilePath);

if ~exist('imDataParams','var') && exist('data','var'),  % 2011.09.22 in case variable is named data
  imDataParams = data; clear data;
end
if isfield(data1,'data'), 
  imDataParams = data1.data;
  clear data1; 
elseif isfield(data1,'imDataParams'), 
  imDataParams = data1.imDataParams;
  clear data1; 
end  
% Get handle of default panel content
%  h1 = TabHandles{size(TabHandles,1),1};

%   Delete the previous panel content
%  if ishandle(h1)
%    delete(h1);             % Delete the default content
%  else
try
  delete(hImageAxes);     % Delete the previous image
catch 
  disp('cannot delete')
end
%  end

%   Set the axes and display the image    
ImgSize = round(0.6*PanelHeight);
% $$$ dx = round(PanelWidth/2-ImgSize/2);
% $$$ dy = round(PanelHeight/2-ImgSize/2);;
dx = max(round(PanelWidth/2-ImgSize/2),365);
dy = 30;

hImageAxes = axes('Parent', TabHandles{1,1}, ...
                  'Units', 'pixels', ...
                  'Position', [dx dy ImgSize ImgSize]);

%imDataParams.images = double(imDataParams.images);
imsos = sqrt(sum(abs(imDataParams.images(:,:,ceil(end/2),:,1)),4)); 
sc = max(imsos(:));
imagesc(imsos,'Parent', hImageAxes,[0 0.75*sc]);axis equal tight off;colormap gray; freezeColors;
fs = 24;
title(hImageAxes,'Source image (first TE)','FontSize',fs);


% $$$   %   Make Image Tab active
% $$$   TabSelectCallback(0,0,1);

DisplayData = cell(5,1);

%  imDataParams.images = coilCombine(double(imDataParams.images(:,:,ceil(end/2),:,:)));

[sx,sy,sz,nc,nt] = size(imDataParams.images);

showTEs = '';
for kt=1:nt-1
  showTEs = [showTEs num2str(imDataParams.TE(kt)*1000,'%0.1f') ', '];
end
showTEs = [showTEs num2str(imDataParams.TE(nt)*1000,'%0.1f') ' ms'];



DisplayData{1} = [num2str(sx) ' X ' num2str(sy) ' X ' num2str(sz)];
DisplayData{2} = num2str(nc);
DisplayData{3} = num2str(nt);
DisplayData{4} = showTEs;
DisplayData{5} = [num2str(imDataParams.FieldStrength,'%1.1f') ' T'];
DisplayData{6} = [num2str(imDataParams.PrecessionIsClockwise)];


ColumnNames = {' Dimensions ' ' Coils ' ' Num Echoes ' ' Echo times ' ' Field Strength ' ' Clockwise? '};
Width = 164;
ColumnWidths = {Width Width Width Width Width Width};

%   Create the table
tabwidth = Width*1+169;
tabheight = 112;
dx = 30;%round(PanelWidth-ImgSize-tabwidth - 80);
dy = 30;


uitable('Position',...
        [dx dy tabwidth tabheight],...
        'Parent', TabHandles{1,1}, ...
        'RowName', ColumnNames,...
        'ColumnName', [],...
        'ColumnWidth', ColumnWidths,...
        'Data', DisplayData);


% $$$   cbh = uicontrol(fh,'Style','checkbox',...
% $$$                 'String','Coil combine before processing',...
% $$$                 'Value',1,'Position',[dx+tabwidth dy+40 130 20]);


% Enable processing button, disable refinement and postprocess
set(butProcess,'enable','on');
set(butProcessRef,'enable','off');
set(butPostProcess,'enable','off');
set(butSave,'enable','off');
set(butSaveRef,'enable','off');
set(butSavePost,'enable','off');


% Clear previous results and images, if any
clear outParams outParamsRef
try 
  delete(wfAxes);
end
try 
  delete(wfAxesRef);
end
try 
  delete(postAxes);
end


%%%% STOP WAIT %%%%%
stopWait( 1 );
stopWait( 2 );
stopWait( 3 );
stopWait( 4 );
%%%%%%%%%%%%%%%%%%%%%



end



function setAlgo(hObj,event)

TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};


global allParams algoParams tab_algoParams1 tab_algoParams2 algoName butProcess algoRowNames ndb mytext2;

% Called when user activates popup menu 
val = get(hObj,'Value');

algoParams = allParams(val).algoParams;
algoName = allParams(val).algoName;

DisplayData = cell(6,1);

fstring1 = [''];
rstring1 = [''];
for kn=1:length(algoParams.species(1).frequency)
  fstring1 = [fstring1 num2str(algoParams.species(1).frequency(kn),'%0.3f') ' '];
  rstring1 = [rstring1 num2str(algoParams.species(1).relAmps(kn),'%0.3f') ' '];
end

fstring2 = [''];  
rstring2 = [''];  
for kn=1:length(algoParams.species(2).frequency)
  fstring2 = [fstring2 num2str(algoParams.species(2).frequency(kn),'%0.3f') ' '];
  rstring2 = [rstring2 num2str(algoParams.species(2).relAmps(kn),'%0.3f') ' '];
end

DisplayData{1} = algoParams.species(1).name;
DisplayData{2} = fstring1;
DisplayData{3} = rstring1;
DisplayData{4} = algoParams.species(2).name;
DisplayData{5} = fstring2;
DisplayData{6} = rstring2;

ColumnNames = {' Species 1 ' ' Frequency ' ' Relative amplitudes  ' ' Species 2 ' ' Frequency ' ' Relative amplitudes  '};
Width = 340;
ColumnWidths = {Width Width Width Width Width Width};
ColEdit = [true true true true true true];

%   Create the table
tabwidth = Width*1+224;
tabheight = 130;
dx = 20;
dy = 80;

try
  delete(tab_algoParams1);
end
tab_algoParams1 = uitable('Position',...
                          [dx dy tabwidth tabheight],...
                          'Parent', TabHandles{2,1}, ...
                          'RowName', ColumnNames,...
                          'ColumnName', [],...
                          'ColumnWidth', ColumnWidths,...
                          'CellEditCallback', @speciesEditCallback,...
                          'ColumnEditable', ColEdit,...
                          'Data', DisplayData);  


Width = 140;
%ColumnWidths = {Width Width Width Width Width Width};

clear DisplayData ColumnWidths;
names = fieldnames(algoParams);
num_rows = length(names)-1;


algoRowNames = [];
ColumnWidths = [];
DisplayData = [];
try
  delete(tab_algoParams2);
end
tabheight=0;
if num_rows>0
  DisplayData = cell(num_rows,1);
  k2 = 0;
  clear DisplayData;
  clear ColEdit;
  for k=1:length(names)
    if ~strcmp(names(k),'species')
      k2 = k2+1;
      algoRowNames{k2} = names{k};
      DisplayData{k2,1} = num2str(getfield(algoParams,names{k}));
      ColumnWidths{k2} = Width;
      ColEdit(k2) = true;
    end
  end
  
  
  tabwidth = Width*1+224;
  tabheight = 21*num_rows + 1;
  dx = 20;
  dy = 220;
  
  tab_algoParams2 = uitable('Position',...
                            [dx dy tabwidth tabheight],...
                            'Parent', TabHandles{2,1}, ...
                            'RowName', algoRowNames,...
                            'ColumnName', [],...
                            'ColumnEditable', ColEdit,...
                            'CellEditCallback', @algoEditCallback,...
                            'ColumnWidth', ColumnWidths,...
                            'Data', DisplayData);  
  
  
end

try 
  delete(mytext2)
end
mytext2 = uicontrol('Style','text',...
                    'Parent', TabHandles{2,1}, ...
                    'FontSize',14,...
                    'String','View/edit parameters',...
                    'Position',[20 225+tabheight tabwidth 40]);




% Button for fat spectrum generation
%%   Define content for the Open Image File Tab
%   Open Image Pushbutton
butModel = uicontrol('Parent', TabHandles{2,1}, ...
                     'Units', 'pixels', ...
                     'Position', [20 20 160 40], ...
                     'String', 'Get fat spectrum', ...
                     'Callback', @getFatSpectrum , ...
                     'Style', 'pushbutton',...
                     'HorizontalAlignment', 'center',...
                     'FontName', 'arial',...
                     'Enable', 'on',...
                     'FontWeight', 'bold',...
                     'FontSize', 12);

ndb = 2;
ColWidths{1} = 100;
ColWidths{2} = 100;
DisplayNDB{1} = num2str(ndb);
ColEditNDB(1) = true;
rowNameNDB{1} = 'Number of double bonds';

tab_ndb = uitable('Position',...
                  [190 20 332 40],...
                  'Parent', TabHandles{2,1}, ...
                  'RowName', rowNameNDB,...
                  'ColumnName', [],...
                  'ColumnEditable', ColEditNDB,...
                  'CellEditCallback', @fatNDBCallback,...
                  'ColumnWidth', ColWidths,...
                  'Data', DisplayNDB);    



%%%% STOP WAIT %%%%%
stopWait( 2 );
%%%%%%%%%%%%%%%%%%%%%


end



function setAlgoRef(hObj,event)

TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};


global allRefParams algoParamsRef tab_algoParamsRef1 tab_algoParamsRef2 algoNameRef butProcess algoRowNames ndb mytext3;


% Called when user activates popup menu 
val = get(hObj,'Value');

algoParamsRef = allRefParams(val).algoParams;
algoNameRef = allRefParams(val).algoName;



DisplayData = cell(6,1);

fstring1 = [''];
rstring1 = [''];
for kn=1:length(algoParamsRef.species(1).frequency)
  fstring1 = [fstring1 num2str(algoParamsRef.species(1).frequency(kn),'%0.3f') ' '];
  rstring1 = [rstring1 num2str(algoParamsRef.species(1).relAmps(kn),'%0.3f') ' '];
end

fstring2 = [''];  
rstring2 = [''];  
for kn=1:length(algoParamsRef.species(2).frequency)
  fstring2 = [fstring2 num2str(algoParamsRef.species(2).frequency(kn),'%0.3f') ' '];
  rstring2 = [rstring2 num2str(algoParamsRef.species(2).relAmps(kn),'%0.3f') ' '];
end

DisplayData{1} = algoParamsRef.species(1).name;
DisplayData{2} = fstring1;
DisplayData{3} = rstring1;
DisplayData{4} = algoParamsRef.species(2).name;
DisplayData{5} = fstring2;
DisplayData{6} = rstring2;

ColumnNames = {' Species 1 ' ' Frequency ' ' Relative amplitudes  ' ' Species 2 ' ' Frequency ' ' Relative amplitudes  '};
Width = 340;
ColumnWidths = {Width Width Width Width Width Width};
ColEdit = [true true true true true true];

%   Create the table
tabwidth = Width*1+224;
tabheight = 130;
dx = 20;
dy = 80;

try
  delete(tab_algoParamsRef1);
end
tab_algoParamsRef1 = uitable('Position',...
                          [dx dy tabwidth tabheight],...
                          'Parent', TabHandles{3,1}, ...
                          'RowName', ColumnNames,...
                          'ColumnName', [],...
                          'ColumnWidth', ColumnWidths,...
                          'CellEditCallback', @speciesEditCallback,...
                          'ColumnEditable', ColEdit,...
                          'Data', DisplayData);  


Width = 140;
%ColumnWidths = {Width Width Width Width Width Width};

clear DisplayData ColumnWidths;
names = fieldnames(algoParamsRef);
num_rows = length(names)-1;


algoRowNames = [];
ColumnWidths = [];
DisplayData = [];
try
  delete(tab_algoParamsRef2);
end

tabheight=0;
if num_rows>0
  DisplayData = cell(num_rows,1);
  k2 = 0;
  clear DisplayData;
  clear ColEdit;
  for k=1:length(names)
    if ~strcmp(names(k),'species')
      k2 = k2+1;
      algoRowNames{k2} = names{k};
      DisplayData{k2,1} = num2str(getfield(algoParamsRef,names{k}));
      ColumnWidths{k2} = Width;
      ColEdit(k2) = true;
    end
  end
  
  
  tabwidth = Width*1+224;
  tabheight = 21*num_rows + 1;
  dx = 20;
  dy = 220;
  
  tab_algoParamsRef2 = uitable('Position',...
                               [dx dy tabwidth tabheight],...
                               'Parent', TabHandles{3,1}, ...
                               'RowName', algoRowNames,...
                               'ColumnName', [],...
                               'ColumnEditable', ColEdit,...
                               'CellEditCallback', @algoEditCallback,...
                               'ColumnWidth', ColumnWidths,...
                               'Data', DisplayData);  
  
  
  
%%%% STOP WAIT %%%%%
stopWait( 3 );
%%%%%%%%%%%%%%%%%%%%%

end

try 
  delete(mytext3)
end
mytext3 = uicontrol('Style','text',...
                    'Parent', TabHandles{3,1}, ...
                    'FontSize',14,...
                    'String','View/edit parameters',...
                    'Position',[20 225+tabheight tabwidth 40]);

% $$$ % Button for fat spectrum generation
% $$$ %%   Define content for the Open Image File Tab
% $$$ %   Open Image Pushbutton
% $$$ butProcess = uicontrol('Parent', TabHandles{3,1}, ...
% $$$                        'Units', 'pixels', ...
% $$$                        'Position', [20 20 160 40], ...
% $$$                        'String', 'Get fat spectrum', ...
% $$$                        'Callback', @getFatSpectrum , ...
% $$$                        'Style', 'pushbutton',...
% $$$                        'HorizontalAlignment', 'center',...
% $$$                        'FontName', 'arial',...
% $$$                        'Enable', 'on',...
% $$$                        'FontWeight', 'bold',...
% $$$                        'FontSize', 12);
% $$$ 
% $$$ ndb = 2;
% $$$ ColWidths{1} = 100;
% $$$ ColWidths{2} = 100;
% $$$ DisplayNDB{1} = num2str(ndb);
% $$$ ColEditNDB(1) = true;
% $$$ rowNameNDB{1} = 'Number of double bonds';
% $$$ 
% $$$ tab_ndb = uitable('Position',...
% $$$                   [190 28 332 25],...
% $$$                   'Parent', TabHandles{3,1}, ...
% $$$                   'RowName', rowNameNDB,...
% $$$                   'ColumnName', [],...
% $$$                   'ColumnEditable', ColEditNDB,...
% $$$                   'CellEditCallback', @fatNDBCallback,...
% $$$                   'ColumnWidth', ColWidths,...
% $$$                   'Data', DisplayNDB);    

end



function setAlgoPost(hObj,event)

TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};


global allPostParams postParams tab_postParams1 algoPostRowNames mytext4;

% Called when user activates popup menu 
val = get(hObj,'Value');
postParams = allPostParams(val).postParams;


clear DisplayData ColumnWidths;
names = fieldnames(postParams);
num_rows = length(names);

Width = 140;
tabwidth = Width*1+224;
tabheight = 22*num_rows + 15;
dx = 20;
dy = 210;

algoRowNames = [];
ColumnWidths = [];
DisplayData = [];
try
  delete(tab_postParams1);
end


if num_rows>0
  DisplayData = cell(num_rows,1);
  k2 = 0;
  clear DisplayData;
  clear ColEdit;
  for k=1:length(names)
    k2 = k2+1;
    algoPostRowNames{k2} = names{k};
    DisplayData{k2,1} = num2str(getfield(postParams,names{k}));
    ColumnWidths{k2} = Width;
    ColEdit(k2) = true;
  end
  
  
  
  tab_postParams1 = uitable('Position',...
                            [dx dy tabwidth tabheight],...
                            'Parent', TabHandles{4,1}, ...
                            'RowName', algoPostRowNames,...
                            'ColumnName', [],...
                            'ColumnEditable', ColEdit,...
                            'CellEditCallback', @algoPostEditCallback,...
                            'ColumnWidth', ColumnWidths,...
                            'Data', DisplayData);  
  

  
try 
  delete(mytext4)
end
mytext4 = uicontrol('Style','text',...
                    'Parent', TabHandles{4,1}, ...
                    'FontSize',14,...
                    'String','View/edit parameters',...
                    'Position',[20 230+tabheight tabwidth 40]);
  
  
end



%%%% STOP WAIT %%%%%
stopWait( 4 );
%%%%%%%%%%%%%%%%%%%%%

end



function getFatSpectrum(ev,source)


TabHandles = guidata(gcf);
NumberOfTabs = size(TabHandles,1)-2;
PanelWidth = TabHandles{NumberOfTabs+1,2};
PanelHeight = TabHandles{NumberOfTabs+1,3};

global ndb algoRowNames algoParams

waterFrequencyPPM = 4.7;
algoParams.species = getFatWaterSpectrum_massScale_fatByDoubleBonds(ndb,waterFrequencyPPM);

global allParams algoParams tab_algoParams1 tab_algoParams2 algoName butProcess algoRowNames;

DisplayData = cell(6,1);

fstring1 = [''];
rstring1 = [''];
for kn=1:length(algoParams.species(1).frequency)
  fstring1 = [fstring1 num2str(algoParams.species(1).frequency(kn),'%0.3f') ' '];
  rstring1 = [rstring1 num2str(algoParams.species(1).relAmps(kn),'%0.3f') ' '];
end

fstring2 = [''];  
rstring2 = [''];  
for kn=1:length(algoParams.species(2).frequency)
  fstring2 = [fstring2 num2str(algoParams.species(2).frequency(kn),'%0.3f') ' '];
  rstring2 = [rstring2 num2str(algoParams.species(2).relAmps(kn),'%0.3f') ' '];
end

DisplayData{1} = algoParams.species(1).name;
DisplayData{2} = fstring1;
DisplayData{3} = rstring1;
DisplayData{4} = algoParams.species(2).name;
DisplayData{5} = fstring2;
DisplayData{6} = rstring2;

ColumnNames = {' Species 1 ' ' Frequency ' ' Relative amplitudes  ' ' Species 2 ' ' Frequency ' ' Relative amplitudes  '};
Width = 340;
ColumnWidths = {Width Width Width Width Width Width};
ColEdit = [true true true true true true];

%   Create the table
tabwidth = Width*1+224;
tabheight = 130;
dx = 20;
dy = 80;

try
  delete(tab_algoParams1);
end
tab_algoParams1 = uitable('Position',...
                          [dx dy tabwidth tabheight],...
                          'Parent', TabHandles{2,1}, ...
                          'RowName', ColumnNames,...
                          'ColumnName', [],...
                          'ColumnWidth', ColumnWidths,...
                          'CellEditCallback', @speciesEditCallback,...
                          'ColumnEditable', ColEdit,...
                          'Data', DisplayData);  


end


function fatNDBCallback(source,ev)

global ndb

currow = ev.Indices(1);

olddata = ev.PreviousData;
newdata = ev.NewData;

newval = str2num(newdata);
if ~isempty(newval)
  ndb = newval;
end

end



function algoEditCallback(source,ev)

global algoRowNames algoParams

currow = ev.Indices(1);

olddata = ev.PreviousData;
newdata = ev.NewData;

newval = str2num(newdata);
if ~isempty(newval)
  algoParams = setfield(algoParams,algoRowNames{currow},newval);
end

end


function algoPostEditCallback(source,ev)

global algoPostRowNames postParams

currow = ev.Indices(1);

olddata = ev.PreviousData;
newdata = ev.NewData;

newval = str2num(newdata);
if ~isempty(newval)
  postParams = setfield(postParams,algoPostRowNames{currow},newval);
end

end


function speciesEditCallback(source,ev)


global algoParams 

currow = ev.Indices(1);

olddata = ev.PreviousData;
newdata = ev.NewData;


switch currow
 case 1
  algoParams.species(1).name = newdata;
 case 2
  newval = str2num(newdata);
  if ~isempty(newval)
    algoParams.species(1).frequency = newval;
  end  
 case 3
  newval = str2num(newdata);
  if ~isempty(newval)
    algoParams.species(1).relAmps = newval;
  end  
  
 case 4
  algoParams.species(2).name = newdata;
 case 5
  newval = str2num(newdata);
  if ~isempty(newval)
    algoParams.species(2).frequency = newval;
  end  
 case 6
  newval = str2num(newdata);
  if ~isempty(newval)
    algoParams.species(2).relAmps = newval;
  end  
  
end


% $$$ fieldName = data{row, 1};
% $$$ fieldValue  = data{row, 2};
% $$$ 
% $$$ str2num(fieldValue).'




end


function initGlobal( )

global allParams;
global allRefParams;
global imDataParams;
global outParams;
global allPostParams
global algoParams algoName algoParamsRef algoNameRef postParams

% Berglund, 2 point
clear algoParams;
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
algoParams.species(2).relAmps = [88 642 58 62 58 6 39 10 37];

allParams(1).algoParams = algoParams;clear algoParams;
allParams(1).algoName = 'fw_i3cm0c_2flexiblepoint_berglund';

% Berglund, 3 point
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
algoParams.species(2).relAmps = [88 642 58 62 58 6 39 10 37];

% Algorithm-specific parameters
algoParams.c1 = 0.75; % Magnitude weight threshold for seed points
algoParams.c2 = 0.25; % Threshold on |log(W/F)| for seed points
allParams(2).algoParams = algoParams;clear algoParams;
allParams(2).algoName = 'fw_i3cm0i_3point_berglund';

% Tsao/Jiang, 3+ poing, multipeak
algoParams.species(1).name = 'Water';
algoParams.species(1).frequency = 4.7;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'Fat (6 peaks)';
algoParams.species(2).frequency = [0.9000    1.3000    2.1000    2.7600    4.3100    5.3000];  % ppm
algoParams.species(2).relAmps =   [0.0871    0.6937    0.1281    0.0040    0.0390    0.0480];
% $$$ algoParams.species(4).name = 'Acetone';
% $$$ algoParams.species(4).frequency = [1.9000];  % ppm
% $$$ algoParams.species(4).relAmps =   [1.0000];
algoParams.MinFractSizeToDivide = 0.01;
algoParams.MaxNumDiv = 7;
allParams(3).algoParams = algoParams;clear algoParams;
allParams(3).algoName = 'fw_i2cm0c_3pluspoint_tsaojiang';

% Tsao/Jiang, 3 poing, singlepeak
algoParams.species(1).name = 'Water';
algoParams.species(1).frequency = 4.7;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'Fat (1 peak)';
algoParams.species(2).frequency = [1.3000];  % ppm
algoParams.species(2).relAmps =   [1.0000];
% $$$ algoParams.species(3).name = 'Acetone';
% $$$ algoParams.species(3).frequency = [1.9000];  % ppm
% $$$ algoParams.species(3).relAmps =   [1.0000];
algoParams.MinFractSizeToDivide = 0.01;
algoParams.MaxNumDiv = 7;
allParams(4).algoParams = algoParams;clear algoParams;
allParams(4).algoName = 'fw_i2cs0c_3point_tsaojiang';

% Hernando, 3+ point, multipeak
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 0]; % Range of R2* values
algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
algoParams.range_fm = [-700 700]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 0; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.01; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
allParams(5).algoParams = algoParams;clear algoParams;
allParams(5).algoName = 'fw_i2cm1i_3pluspoint_hernando_graphcut';


% Lu, 3 point, multipeak
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
allParams(6).algoParams = algoParams;clear algoParams;
allParams(6).algoName = 'fw_3point_wm_goldSect';


% RG VARIABLES - see Yu, et al. MRM 2005 54:1032-1039 for details
algoParams.reg_grow=1;    % 0 - No (voxel independent), 1 - Yes (RG)
algoParams.rg=[15 15];    % entry 1: number of super pixels               
                          % entry 2: extrapolation / region growing kernel size ...
                          % weighting neighborhood

algoParams.downsize=[16 16];% downsample size of low res image to start RG
algoParams.filtsize=[7 7];  % median filter size for final smoothing of fieldmap
algoParams.maxiter=300;     % maximum number of IDEAL iterations before loop break
algoParams.fm_epsilon=1;    % in IDEAL iteration, the fm error tolerance (Hz)
algoParams.sp_mp=1;         % single fat peak (0) or multi fat peak (1)
                            %% *************************************************************************
                            % SPECTRAL MODEL
                            %% *************************************************************************
algoParams.M = 2;% Number of species; 2 = fat/water
                 % water
algoParams.species(1).name = 'water';
algoParams.species(1).frequency=0; 
algoParams.species(1).relAmps=1;

% fat
algoParams.species(2).name = 'fat';

if algoParams.sp_mp==1    
  % NOTE THIS IS THE DIFFERENCE IN PPM BETWEEN FAT AND WATER, NOT THE
  % ACTUAL FAT PPMs
  
  algoParams.species(2).frequency=-[3.8 3.4 2.6 1.94 0.39 -0.6]; % positive fat frequency --> -2*pi*fieldmap*te effect
  algoParams.species(2).relAmps=[0.087 0.693 0.128 0.004 0.039 0.048];
else
  algoParams.species(2).frequency=-3.4;
  algoParams.species(2).relAmps=1;
end
allParams(7).algoParams = algoParams;clear algoParams;
allParams(7).algoName = 'fw_i2cm0c_3pluspoint_RGsinglecoil_hu';


% Sharma, Multipeak, 3+ point
%* BEGIN: Set algorithm parameters *
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% algorithm specific parameters
algoParams.stepsize = 0.75;     % support size scaling 
algoParams.min_win_size = 16;   % minimum 1D support size for B-spline
algoParams.MaxIter = 10;        % maximum iterations per scale
                                %* END: Set algorithm parameters *
allParams(8).algoParams = algoParams;clear algoParams;
allParams(8).algoName = 'fw_i2cm0i_3pluspoint_sharma';

algoParams = allParams(1).algoParams;
algoName = allParams(1).algoName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Refinements    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Hernando, 3+ point, multipeak, mixed fitting
% General parameters
clear algoParamsRef;
algoParamsRef.species(1).name = 'water';
algoParamsRef.species(1).frequency = 0;
algoParamsRef.species(1).relAmps = 1;
algoParamsRef.species(2).name = 'fat';
algoParamsRef.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParamsRef.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
% Algorithm-specific parameters
algoParamsRef.NUM_MAGN = 1; % Number of echoes to discard phase at beginning of train
algoParamsRef.THRESHOLD = 0.04; % Threshold for processing each voxel (relative to max signal over all voxels)
algoParamsRef.range_r2star = [0 200]; % Range of admissible R2* values
allRefParams(1).algoParams = algoParamsRef;clear algoParamsRef;
allRefParams(1).algoName = 'fw_i2xm1c_3pluspoint_hernando_mixedfit';


% Eggers: spiral deblurring
% Set parameters for field inhomogeneity correction
algoParamsRef.species(1).name = 'water';
algoParamsRef.species(1).frequency = 0;
algoParamsRef.species(1).relAmps = 1;
algoParamsRef.species(2).name = 'fat';
algoParamsRef.species(2).frequency = [-3.4];
algoParamsRef.species(2).relAmps = [1];
algoParamsRef.nv      =  15;         % Number of spiral interleaves
algoParamsRef.ns      =  6434;       % Number of spiral samples per interleaf
algoParamsRef.alpha   =  1.25;       % Gridding oversampling factor
algoParamsRef.ks      =  4;          % Gridding kernel size
algoParamsRef.gti     =  1000;       % Gridding kernel table increment
algoParamsRef.tau     =  2.325e-6;   % Sampling interval [s] 
algoParamsRef.shutter =  1;          % Circular image space shutter
algoParamsRef.lambda  =  3.0;        % Shape parameter
algoParamsRef.delay   =  0.0;
allRefParams(2).algoParams = algoParamsRef;clear algoParamsRef;
allRefParams(2).algoName = 'cpr';

algoParamsRef = allRefParams(1).algoParams;
algoNameRef = allRefParams(1).algoName;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear postParams;
postParams.include_r2star = 0;
allPostParams(1).postParams = postParams; clear postParams;

postParams.noise_bias_correction = 1;
allPostParams(2).postParams = postParams; clear postParams;

postParams = allPostParams(1).postParams;



end
