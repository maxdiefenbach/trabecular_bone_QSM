function varargout = fw_inputparams(varargin)
% FW_INPUTPARAMS M-file for fw_inputparams.fig
%      FW_INPUTPARAMS, by itself, creates a new FW_INPUTPARAMS or raises the existing
%      singleton*.
%
%      H = FW_INPUTPARAMS returns the handle to a new FW_INPUTPARAMS or the handle to
%      the existing singleton*.
%
%      FW_INPUTPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FW_INPUTPARAMS.M with the given input arguments.
%
%      FW_INPUTPARAMS('Property','Value',...) creates a new FW_INPUTPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fw_inputparams_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fw_inputparams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fw_inputparams

% Last Modified by GUIDE v2.5 10-Nov-2011 18:57:34

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ... %YunJiang20111221
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fw_inputparams_OpeningFcn, ...
                   'gui_OutputFcn',  @fw_inputparams_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function species = DefaultLibrary()
  % Default
  species(1).name = 'Water';
  species(1).frequency = 4.7;
  species(1).relAmps = 1;
  species(2).name = 'Fat (6 peaks)';
  species(2).frequency = [0.9000    1.3000    2.1000    2.7600    4.3100    5.3000];  % ppm
  species(2).relAmps =   [0.0871    0.6937    0.1281    0.0040    0.0390    0.0480];
  species(3).name = 'Fat (1 peak)';
  species(3).frequency = [1.3000];  % ppm
  species(3).relAmps =   [1.0000];
  species(4).name = 'Acetone';
  species(4).frequency = [1.9000];  % ppm
  species(4).relAmps =   [1.0000];


function  status = WriteLibrary(LibraryFilename,species)
  % Yun Jiang - Debug
  if 0, fprintf('In function of WriteLibrary\n');end;
  [tmpdir,tmpfile] = fileparts(LibraryFilename);
  if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
  if ~exist(tmpdir,'dir'), status = 1; return; end % Error
  clear tmpdir tmpfile;
  
  fid=-1;
  try,
    fid = fopen(LibraryFilename,'w');
    if fid<0, status = 2; return; end % Error
    for n=1:numel(species),
      if ~isempty(species(n).name)
        fprintf(fid,'%s\n',species(n).name);
        fprintf(fid,'%f',species(n).frequency(1)); fprintf(fid,'\t%f',species(n).frequency(2:end)); fprintf(fid,'\n');
        fprintf(fid,'%f',species(n).relAmps(1)); fprintf(fid,'\t%f',species(n).relAmps(2:end)); fprintf(fid,'\n');
        fprintf(fid,'\n');
      end
    end; clear n;
    fclose(fid);
    clear fid;
    status=0;
    fprintf('Saved library to %s\n',LibraryFilename);
  catch,
    if fid, try, fclose(fid); end; end;
    clear fid;
    status = 3;
  end

function  species = ReadLibrary(LibraryFilename)
  % Yun Jiang - Debug
  if 0, fprintf('In function of ReadLibrary\n');end;
    species=[]; 
  if ~exist(LibraryFilename,'file'), return; end  % Error
  fid = -1;
  try,
    fid = fopen(LibraryFilename,'r');
    if fid<0, return; end % Error
    NumSpecies = 0;
    while 1
      tmpline = fgetl(fid);
      if ~ischar(tmpline), break; end
      tmpname = tmpline;
      tmpline = fgetl(fid);
      if ~ischar(tmpline), break; end
      tmpfrequency = sscanf(tmpline,'%f');
      tmpline = fgetl(fid);
      if ~ischar(tmpline), break; end
      tmprelAmps = sscanf(tmpline,'%f');
      tmpline = fgetl(fid);
      if ~ischar(tmpline), break; end
      clear tmpline;      
      if numel(tmpfrequency)~=numel(tmprelAmps),
        n = min(numel(tmpfrequency),numel(tmprelAmps));
        tmpfrequency = tmpfrequency(1:n);
        tmprelAmps = tmprelAmps(1:n);
        clear n;
      end
      if ~isempty(tmpname) & numel(tmpfrequency)>=1
        NumSpecies = NumSpecies+1;
        species(NumSpecies).name = tmpname;
        species(NumSpecies).frequency = tmpfrequency;
        species(NumSpecies).relAmps = tmprelAmps;
      end
    end
    clear tmpname tmpfrequency tmpAmps NumSpecies;
    fclose(fid);
    clear fid;
  catch,
    if fid>0, try, fclose(fid); end; end; clear fid;
    return; 
  end
  
function handles = PopulateLibraryList(hObject, eventdata, handles, species)
  % Yun Jiang - Debug
  if 0, fprintf('In function of PopulateLibraryList\n');end;
    NumSpecies = 0;
  for n=1:numel(species),
    if ~isempty(species(n).name) & numel(species(n).frequency)>=1 & numel(species(n).frequency)==numel(species(n).relAmps),
      NumSpecies = NumSpecies + 1;
    end
  end; clear n;
  handles.species = [];
  handles.library_toppos = 1;
  handles.library_selection = [1,2];
  NumSpecies = 0;
  for n=1:numel(species),
    if ~isempty(species(n).name) & numel(species(n).frequency)>=1 & numel(species(n).frequency)==numel(species(n).relAmps),
      NumSpecies = NumSpecies + 1;
      handles.species(NumSpecies).name = species(n).name;
      handles.species(NumSpecies).frequency = species(n).frequency;
      handles.species(NumSpecies).relAmps = species(n).relAmps;
    end
  end; clear n;
  %if evalin('caller','exist(''handles'',''var'')') assignin('caller','handles',handles); end
  %guidata(hObject, handles); % Update handles structure  
  UpdateLibraryList(hObject, eventdata, handles);

function UpdateLibraryList(hObject, eventdata, handles)
    % Yun Jiang - Debug
  if 0, fprintf('In function of UpdateLibraryList\n');end;
  if numel(handles.species)+1 <= numel(handles.checkboxes),
    set(handles.slider1,'enable','off','visible','off');
  else
    maxval = numel(handles.species)+1-numel(handles.checkboxes)+1;
    set(handles.slider1,'Min',1,'Max',maxval,'Value',maxval-handles.library_toppos+1,...
      'enable','on','visible','on',...
      'sliderstep',[1/(maxval-1),5/(maxval-1)]);
    clear maxval;
  end
  for n=1:numel(handles.checkboxes),
    NumInList = handles.library_toppos+n-1;
    if NumInList <= numel(handles.species)
      set(handles.checkboxes(n),'String',handles.species(NumInList).name,'enable','on','visible','on','FontAngle','normal',...
        'Value',any(handles.library_selection == NumInList));
    elseif NumInList == numel(handles.species)+1
      set(handles.checkboxes(n),'String','+ New species...','enable','on','visible','on','FontAngle','italic','Value',0);
    else
      set(handles.checkboxes(n),'enable','off','visible','off');
    end
    clear NumInList;
  end; clear n;
  
% --- Executes just before fw_inputparams is made visible.
function fw_inputparams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fw_inputparams (see VARARGIN)

% Yun Jiang - Debug
  if 0, fprintf('In function of fw_inputparams_OpeningFcn\n');end;

% Add path 
[BASEPATH,tmpfile] = fileparts(mfilename('fullpath'));
tmp = BASEPATH; addpath(tmp); fprintf('Adding to path: %s\n',tmp); clear tmp;
handles.LibraryFilename = fullfile(BASEPATH,'private',[tmpfile,'.txt']); clear tmpfile;

% GUI
handles.checkboxes = [handles.checkbox1,handles.checkbox2,handles.checkbox3,handles.checkbox4,handles.checkbox5,...
                      handles.checkbox6,handles.checkbox7,handles.checkbox8,handles.checkbox9,handles.checkbox10];
handles.speciesdetails = [handles.text_speciesname,handles.text_ppm,handles.text_relAmps,...
                          handles.edit_speciesname,handles.edit_ppm,handles.edit_relAmps,...
                          handles.text_status,...
                          handles.pushbutton_deletespecies,...
                          handles.pushbutton_savespecies,...
                          handles.pushbutton_cancelspecies];
                        % Reset
handles.warningstatus = [handles.uipanel_library, ...
                         handles.text_speciesname, ...
                         handles.text_ppm, ...
                         handles.text_relAmps, ...
                         handles.text_fieldstrength, ...
                         handles.text_precession];
  
handles.speciesdetails_current = [];
handles.speciesdetails_changed = 0;
for n=1:numel(handles.speciesdetails),
  set(handles.speciesdetails(n),'visible','off');
end; clear n;

% Library
if ~exist(handles.LibraryFilename,'file'),
  species = [];
else
  species = ReadLibrary(handles.LibraryFilename);
end
if isempty(species)
  species = DefaultLibrary();
  WriteLibrary(handles.LibraryFilename,species);
  fprintf('Spectral library created: %s\n',handles.LibraryFilename);
end
handles = PopulateLibraryList(hObject, eventdata, handles, species);
tmpselection = handles.library_selection;
handles.library_selection = [];

for arginnum=1:numel(varargin)
  if isfield(varargin{arginnum},'species'),  % Species
    for n=1:numel(varargin{arginnum}.species),
      if isfield(varargin{arginnum}.species(n),'name'),
        tmpname = strtrim(varargin{arginnum}.species(n).name);
        if isfield(varargin{arginnum}.species(n),'frequency'),
          tmpfrequency = varargin{arginnum}.species(n).frequency;
          if isfield(varargin{arginnum}.species(n),'relAmps'),
            tmprelAmps = varargin{arginnum}.species(n).relAmps;
            Found = 0;
            for m=1:numel(handles.species)
              if isequal(handles.species(m).name     ,varargin{arginnum}.species(n).name) ...
               & isequal(handles.species(m).frequency,varargin{arginnum}.species(n).frequency) ...
               & isequal(handles.species(m).relAmps  ,varargin{arginnum}.species(n).relAmps),
                Found = 1;% found
                break;
              end
            end
            if Found, 
              handles.library_selection(end+1) = m;
            else
              m = numel(handles.species)+1;
              handles.species(m).name      = varargin{arginnum}.species(n).name;
              handles.species(m).frequency = varargin{arginnum}.species(n).frequency;
              handles.species(m).relAmps   = varargin{arginnum}.species(n).relAmps;
              handles.library_selection(end+1) = m;
            end
            clear m Found tmprelAmps;
          end; clear tmpfrequency;
        end; clear tmpname;
      end
    end; clear n;
  end
  if isfield(varargin{arginnum},'PrecessionIsClockwise'),  % Species
    tmpvalue = varargin{arginnum}.PrecessionIsClockwise;
    if isempty(tmpvalue) | ~isnumeric(tmpvalue) | tmpvalue>0,
      % Clockwise
      tmpcontents = get(handles.popupmenu_precession,'String');
      for n=1:numel(tmpcontents)
        if isequal(tmpcontents{n},'Clockwise'),
          set(handles.popupmenu_precession,'Value',n);
          break;
        end
      end; clear n tmpcontents;
    else
      % Counter clockwise
      tmpcontents = get(handles.popupmenu_precession,'String');
      for n=1:numel(tmpcontents)
        if ~isequal(tmpcontents{n},'Clockwise'),
          set(handles.popupmenu_precession,'Value',n);
          break;
        end
      end; clear n tmpcontents;
    end; clear tmpvalue;
  end
  if isfield(varargin{arginnum},'FieldStrength'),  % Species
    tmpvalue = varargin{arginnum}.FieldStrength;
    if ~isempty(tmpvalue) & isnumeric(tmpvalue) & isreal(tmpvalue) & tmpvalue>=0 & tmpvalue<=1000,  % Valid
      tmpstr = sprintf('%g',tmpvalue(1));
      set(handles.edit_fieldstrength, 'String', tmpstr, 'UserData', tmpstr);
    end; clear tmpvalue;
  end
end; clear arginnum;

handles.library_selection = unique(handles.library_selection);
if isempty(handles.library_selection), handles.library_selection=[1,2]; end
if ~isequal(tmpselection,handles.library_selection),
  UpdateLibraryList(hObject, eventdata, handles);
end;
clear tmpselection;

% Choose default command line output for fw_inputparams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fw_inputparams wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = fw_inputparams_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Yun Jiang - Debug
  if 0, fprintf('In function of fw_inputparams_OutputFcn\n');end;


% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit_ppm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ppm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

% Yun Jiang - Debug
  if 0, fprintf('In function of edit_ppm_CreateFcn\n');end;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_deletespecies.
function pushbutton_deletespecies_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_deletespecies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Yun Jiang - Debug
  if 0, fprintf('In function of pushbutton_deletespecies_Callback\n');end;

  resetwarnings(hObject, eventdata, handles); % Reset warnings

  if numel(handles.species)==1,
    errordlg('Cannot delete last species in library.','Error','modal');
    return;
  end
  if isempty(handles.speciesdetails_current), return; end   % Window not open
  
  tmpstr = strtrim(get(handles.edit_speciesname,'String'));
  if ~isempty(tmpstr), tmpstr = sprintf('''%s''',tmpstr); else tmpstr = 'this species'; end
  status = questdlg({sprintf('Delete %s?',tmpstr),'This CANNOT be undone.'}, 'Confirm',  'Delete','Cancel','Cancel');
  clear tmpstr;
  if isequal(status,'Delete')
    % Remove species and library_selection
    handles.species = handles.species([1:handles.speciesdetails_current-1,handles.speciesdetails_current+1:numel(handles.species)]);
    handles.library_selection = setdiff(handles.library_selection, handles.speciesdetails_current);
    tmpidx = find(handles.library_selection > handles.speciesdetails_current);
    handles.library_selection(tmpidx) = handles.library_selection(tmpidx)-1; clear tmpidx;

    % Update checkbox list
    WhichCheckbox = handles.speciesdetails_current - handles.library_toppos + 1; % handles.speciesdetails_current is in the new entries
    if WhichCheckbox>=1 & WhichCheckbox<=numel(handles.checkboxes),  % Species to be deleted is within current visible checkboxes
      if handles.library_toppos + numel(handles.checkboxes) - 1 > numel(handles.species)+1,
        handles.library_toppos = max(handles.library_toppos-1, 1);
      end
      UpdateLibraryList(hObject, eventdata, handles);
      uicontrol(handles.checkboxes(WhichCheckbox));
    end
    clear WhichCheckbox;

    % Save library
    status = WriteLibrary(handles.LibraryFilename,handles.species);
    if status,
      errordlg({'Cannot save library file',sprintf('%s',handles.LibraryFilename)},'Error','modal');
    end; clear status;
    
    % Close details window
    handles = CloseSpeciesDetailsWindow(hObject, eventdata, handles);
    guidata(hObject, handles); % Update handles structure  
  end
  clear status;


% --- Executes on selection change in popupmenu_precession.
function popupmenu_precession_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_precession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_precession contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_precession

% Yun Jiang - Debug
  if 0, fprintf('In function of popupmenu_precession_Callback\n');end;

    resetwarnings(hObject, eventdata, handles); % Reset warnings


% --- Executes during object creation, after setting all properties.
function popupmenu_precession_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_precession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

% Yun Jiang - Debug
  if 0, fprintf('In function of popupmenu_precession_CreateFcn\n');end;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fieldstrength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fieldstrength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fieldstrength as text
%        str2double(get(hObject,'String')) returns contents of edit_fieldstrength as a double


% Yun Jiang - Debug
  if 0, fprintf('In function of edit_fieldstrength_Callback\n');end;

  resetwarnings(hObject, eventdata, handles); % Reset warnings

  tmpvalue = str2double(get(hObject,'String'));
  if isempty(tmpvalue) | isnan(tmpvalue),
    set(hObject,'String',get(hObject,'UserData'));
  elseif tmpvalue<0 | tmpvalue>1000
    set(handles.text_fieldstrength,'Foreground',[1,0,0],'FontWeight','bold');
    set(hObject,'String',get(hObject,'UserData'));
  else
    tmpstr = sprintf('%g',tmpvalue(1));
    set(hObject,'String',tmpstr,'UserData',tmpstr);
    clear tmpstr;
  end
  clear tmpvalue;
    

% --- Executes during object creation, after setting all properties.
function edit_fieldstrength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fieldstrength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

% Yun Jiang - Debug
  if 0, fprintf('In function of edit_fieldstrength_CreateFcn\n');end;
  
  
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Yun Jiang - Debug
  if 0, fprintf('In function of pushbutton_OK_Callback\n');end;
  
  
  resetwarnings(hObject, eventdata, handles); % Reset warnings
  NotReadyToExit = 0;  
  [handles, Cancel] = FinishSpeciesDetails(hObject,eventdata,handles);
  if Cancel, clear Cancel; return; end; clear Cancel;
  
  if numel(handles.library_selection)<2, 
    NotReadyToExit = 1;
    set(handles.uipanel_library, 'Foreground',[1,0,0],'FontWeight','bold');
  end

  if isempty(get(handles.edit_fieldstrength,'String')),
    NotReadyToExit = 1;
    set(handles.text_fieldstrength, 'Foreground',[1,0,0],'FontWeight','bold');
    uicontrol(handles.edit_fieldstrength); % Re-focus
  end
  
  if NotReadyToExit, clear NotReadyToExit; return; end
  clear NotReadyToExit;
  
  % Collect output parameters
  handles.output = [];
  handles.output.FieldStrength = str2double(get(handles.edit_fieldstrength,'String'));
  for n=1:numel(handles.library_selection),
    handles.output.species(n) = struct('name',      handles.species(handles.library_selection(n)).name, ...
                                       'frequency', handles.species(handles.library_selection(n)).frequency, ...
                                       'relAmps',   handles.species(handles.library_selection(n)).relAmps);
  end; clear n;
  tmpcontents = get(handles.popupmenu_precession,'String');
  switch tmpcontents{get(handles.popupmenu_precession,'Value')}
    case 'Clockwise'
      handles.output.PrecessionIsClockwise = 1;
    case 'Counter-clockwise'
      handles.output.PrecessionIsClockwise = 0;
  end; clear tmpcontents;
  
  guidata(hObject,handles);
  uiresume(handles.figure1); 

% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Yun Jiang - Debug
  if 0, fprintf('In function of pushbutton_Cancel_Callback\n');end;
  
  handles.output = [];
  guidata(hObject,handles);
  uiresume(handles.figure1);
  
% --- Executes during object creation, after setting all properties.
function edit_relAmps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_relAmps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

% Yun Jiang - Debug
  if 0, fprintf('In function of edit_relAmps_CreateFcn\n');end;
  
  
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_speciesname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_speciesname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Yun Jiang - Debug
  if 0, fprintf('In function of edit_speciesname_Callback\n');end;

  resetwarnings(hObject, eventdata, handles); % Reset warnings

  if isequal(get(handles.figure1,'CurrentCharacter'),27),  % escape
    set(hObject,'String',get(hObject,'UserData'));
    return
  end
  
  handles.speciesdetails_changed = 1;
  [CanSave,tmpspeciesname,tmpppm,tmprelAmps] = CanSaveSpeciesDetails(hObject,eventdata,handles);
  guidata(hObject, handles); % Update handles structure

function resetwarnings(hObject, eventdata, handles)
    
    % Yun Jiang - Debug
  if 0, fprintf('In function of resetwarnings\n');end;
  for n=1:numel(handles.warningstatus)
    set(handles.warningstatus(n),'Foreground',[0,0,0],'FontWeight','normal');
  end; clear n;
    
function edit_ppm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_speciesname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 % Yun Jiang - Debug
  if 0, fprintf('In function of edit_ppm_Callback\n');end;

  resetwarnings(hObject, eventdata, handles); % Reset warnings

  if isequal(get(handles.figure1,'CurrentCharacter'),27),  % escape
    set(hObject,'String',get(hObject,'UserData'));
    return
  end
  handles.speciesdetails_changed = 1;
  %set(handles.text_status,'String','');
  %set(handles.text_speciesname,'Foreground',[0,0,0],'FontWeight','normal');
  %set(handles.text_ppm,'Foreground',[0,0,0],'FontWeight','normal');
  %set(handles.text_relAmps,'Foreground',[0,0,0],'FontWeight','normal');
  [CanSave,tmpspeciesname,tmpppm,tmprelAmps] = CanSaveSpeciesDetails(hObject,eventdata,handles);
  guidata(hObject, handles); % Update handles structure

function edit_relAmps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_relAmps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Yun Jiang - Debug
  if 0, fprintf('In function of edit_relAmps_Callback\n');end;

% Hints: get(hObject,'String') returns contents of edit_relAmps as text
%        str2double(get(hObject,'String')) returns contents of edit_relAmps as a double
  if isequal(get(handles.figure1,'CurrentCharacter'),27),  % escape
    set(hObject,'String',get(hObject,'UserData'));
    return
  end
  handles.speciesdetails_changed = 1;
  %set(handles.text_status,'String','');
  %set(handles.text_speciesname,'Foreground',[0,0,0],'FontWeight','normal');
  %set(handles.text_ppm        ,'Foreground',[0,0,0],'FontWeight','normal');
  %set(handles.text_relAmps    ,'Foreground',[0,0,0],'FontWeight','normal');
  [CanSave,tmpspeciesname,tmpppm,tmprelAmps] = CanSaveSpeciesDetails(hObject,eventdata,handles);
  guidata(hObject, handles); % Update handles structure

function [status,tmpstr] = edit_speciesname_Verify(hObject, eventdata, handles)
    % Yun Jiang - Debug
  if 0, fprintf('In function of edit_speciesname_Verify\n');end;
    
  tmpstr = strtrim(get(handles.edit_speciesname,'String'));
  if isempty(tmpstr),
    % No name present
    set(handles.text_speciesname,'Foreground',[1,0,0],'FontWeight','bold');
    status = 1;
  else
    set(handles.text_speciesname,'Foreground',[0,0,0],'FontWeight','normal');
    status = 0;
  end

function [status,tmpvalues] = edit_ppm_Verify(hObject, eventdata, handles)
    
        % Yun Jiang - Debug
  if 0, fprintf('In function of edit_ppm_Verify\n');end;
    
  try,
    tmpstr = get(handles.edit_ppm,'String');
    tmpstr1 = strtrim(tmpstr(1,:));
    for n=2:size(tmpstr,1)
      tmpstr1 = [tmpstr1 '\n' strtrim(tmpstr(n,:))];
    end; clear n;
    tmpvalues = strread(sprintf(tmpstr1),'%f','delimiter',',;\t'); clear tmpstr tmpstr1;
    tmpstr = sprintf('%g\n',tmpvalues);
    set(handles.text_ppm,'Foreground',[0,0,0],'FontWeight','normal');
    set(handles.edit_ppm,'String',tmpstr,'UserData',tmpstr);
    status = 0;  % OK
    return;
  end
  status = 1;  tmpvalues = []; % Failed
  set(handles.text_ppm,'Foreground',[1,0,0],'FontWeight','bold');

function [status,tmpvalues] = edit_relAmps_Verify(hObject, eventdata, handles)
    
           % Yun Jiang - Debug
  if 0, fprintf('In function of edit_relAmps_Verify\n');end;
  
  try,
    tmpstr = get(handles.edit_relAmps,'String');
    tmpstr1 = strtrim(tmpstr(1,:));
    for n=2:size(tmpstr,1)
      tmpstr1 = [tmpstr1 '\n' strtrim(tmpstr(n,:))];
    end; clear n;
    tmpvalues = strread(sprintf(tmpstr1),'%f','delimiter',',;\t'); clear tmpstr tmpstr1;
    tmpstr = sprintf('%g\n',tmpvalues);
    set(handles.text_relAmps,'Foreground',[0,0,0],'FontWeight','normal');
    set(handles.edit_relAmps,'String',tmpstr,'UserData',tmpstr);
    status = 0;  % OK
    return;
  end
  status = 1;  tmpvalues = []; % Failed
  set(handles.text_relAmps,'Foreground',[1,0,0],'FontWeight','bold');
  
function [CanSave,tmpspeciesname,tmpppm,tmprelAmps] = CanSaveSpeciesDetails(hObject,eventdata,handles)
  
           % Yun Jiang - Debug
  if 0, fprintf('In function of CanSaveSpeciesDetails\n');end;
  
  CanSave = 1; tmpspeciesname=[]; tmpppm=[]; tmprelAmps=[];
  set(handles.text_ppm    ,'Foreground',[0,0,0],'FontWeight','normal');
  set(handles.text_relAmps,'Foreground',[0,0,0],'FontWeight','normal');
  [status,tmpspeciesname] = edit_speciesname_Verify(hObject, eventdata, handles); if status, CanSave=0; end
  [status,tmpppm        ] =         edit_ppm_Verify(hObject, eventdata, handles); if status, CanSave=0; end
  [status,tmprelAmps    ] =     edit_relAmps_Verify(hObject, eventdata, handles); if status, CanSave=0; end
  if CanSave,
    if numel(tmpppm)~=numel(tmprelAmps) | numel(tmpppm)==0, % Either mismatch, or no entries
      set(handles.text_ppm    ,'Foreground',[1,0,0],'FontWeight','bold');
      set(handles.text_relAmps,'Foreground',[1,0,0],'FontWeight','bold');
      CanSave = 0;
      set(handles.text_status,'String','Please check item(s) in RED.','ForegroundColor',[1,0,0],'FontWeight','bold');
    else
      set(handles.text_status,'String','');
    end
  else
    set(handles.text_status,'String','Please check item(s) in RED.','ForegroundColor',[1,0,0],'FontWeight','bold');
  end
  clear status;
  
function [status,handles] = SaveSpeciesDetails(hObject, eventdata, handles)
    
           % Yun Jiang - Debug
  if 0, fprintf('In function of SaveSpeciesDetails\n');end;
  
  [CanSave,tmpspeciesname,tmpppm,tmprelAmps] = CanSaveSpeciesDetails(hObject,eventdata,handles);
  if CanSave,
    handles.species(handles.speciesdetails_current).name = tmpspeciesname;
    handles.species(handles.speciesdetails_current).frequency = tmpppm;
    handles.species(handles.speciesdetails_current).relAmps = tmprelAmps;
    WhichCheckbox = handles.speciesdetails_current - handles.library_toppos + 1; % handles.speciesdetails_current is in the new entries
    if WhichCheckbox>=1 & WhichCheckbox<=numel(handles.checkboxes),
      set(handles.checkboxes(WhichCheckbox),'String',tmpspeciesname);
    end
    status = WriteLibrary(handles.LibraryFilename,handles.species);
    if status,
      errordlg({'Cannot save library file',sprintf('%s',handles.LibraryFilename)},'Error','modal');
    else
      set(handles.text_status,'String','Library saved.','ForegroundColor',[0,0,0],'FontWeight','normal');
      handles.speciesdetails_changed = 0;
      if handles.speciesdetails_current == numel(handles.species), % Need to update library list
        UpdateLibraryList(hObject, eventdata, handles);
      end
     end
    %if evalin('caller','exist(''handles'',''var'')') assignin('caller','handles',handles); end
    %guidata(hObject, handles); % Update handles structure
  else
    status = 1;
  end
  clear CanSave tmpspeciesname tmpppm tmprelAmps;


% --- Executes during object creation, after setting all properties.
function edit_speciesname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_speciesname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

           % Yun Jiang - Debug
  if 0, fprintf('In function of edit_speciesname_CreateFcn\n');end;
  
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox1.
function checkbox_Callback(hObject, eventdata, handles, WhichCheckbox)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox1

% Yun Jiang - Debug
  if 0, fprintf('In function of checkbox_Callback\n');end;

  resetwarnings(hObject, eventdata, handles); % Reset warnings
  
  handles.CurrentSpecies = handles.library_toppos + WhichCheckbox-1;
  %if isequal(get(handles.speciesdetails(1),'visible'),'off'), % Detail window not shown
  %  if get(hObject,'Value')==0, % not selected -> select it, and show details
  %    set(hObject,'Value',1);
  %  end
  %end
  set(handles.text_status,'String','');

  [handles, Cancel] = FinishSpeciesDetails(hObject,eventdata,handles);
  if Cancel, clear Cancel; return; end; clear Cancel;
    
  if get(hObject,'Value'), % Checked
    handles.library_selection = unique([handles.library_selection, handles.CurrentSpecies]);
    handles = OpenSpeciesDetailsWindow(hObject, eventdata, handles);
  else                     % Uncheck
    %if ~isequal(handles.CurrentSpecies,handles.speciesdetails_current) & handles.CurrentSpecies~=numel(handles.species)+1,
    %  % Show details first
    %  set(hObject,'Value',1);
    %  handles = OpenSpeciesDetailsWindow(hObject, eventdata, handles);
    %else 
      handles.library_selection = setdiff(handles.library_selection, handles.CurrentSpecies);
      handles = CloseSpeciesDetailsWindow(hObject, eventdata, handles);
    %end
  end
  guidata(hObject, handles); % Update handles structure
  clear currentpos

function [handles, Cancel] = FinishSpeciesDetails(hObject,eventdata,handles)
  
    % Yun Jiang - Debug
  if 0, fprintf('In function of FinishSpeciesDetails\n');end;
    
    Cancel = 0;
  if handles.speciesdetails_changed,   % Details have changed
    set(hObject,'Value',1-get(hObject,'Value'));    % Restore setting for now...
    WhichCheckbox = handles.speciesdetails_current - handles.library_toppos + 1; % handles.CurrentSpecies is in the new entries
    if WhichCheckbox>=1 & WhichCheckbox<=numel(handles.checkboxes),
      uicontrol(handles.checkboxes(WhichCheckbox)); % Restore focus on entry corresponding to details entry
    end
    
    [CanSave,tmpspeciesname,tmpppm,tmprelAmps] = CanSaveSpeciesDetails(hObject,eventdata,handles);
    if ~isempty(tmpspeciesname), tmpstr = sprintf('in ''%s''',tmpspeciesname); else tmpstr = 'on the right'; end
    if CanSave,
      status = questdlg(sprintf('Discard changes %s?',tmpstr), 'Confirm',  'Save changes','Discard changes','Cancel','Cancel');
    else
      status = questdlg(sprintf('Discard changes %s?',tmpstr), 'Confirm',  'Discard changes','Cancel','Cancel');
    end; clear tmpstr;
    if     isequal(status, 'Save changes'),          % Save
      [tmpstatus,handles] = SaveSpeciesDetails(hObject, eventdata, handles);
      clear tmpstatus;
    elseif ~isequal(status, 'Discard changes')       % Cancel
      Cancel = 1;
      return; 
    end
    set(hObject,'Value',1-get(hObject,'Value'));    % Restore new setting
    clear tmpspeciesname tmpppm tmprelAmps CanSave;
  else                                              % No details changed
    status = 'Discard changes';
  end
  
  if isequal(status, 'Discard changes'),       % Discard
    % If details showing new entry, deselect it.
    if handles.speciesdetails_current == numel(handles.species)+1,
      WhichCheckbox = handles.speciesdetails_current - handles.library_toppos + 1; % handles.CurrentSpecies is in the new entries
      if WhichCheckbox>=1 & WhichCheckbox<=numel(handles.checkboxes),
        set(handles.checkboxes(WhichCheckbox),'Value',0);
      end
      handles.library_selection = setdiff(handles.library_selection, handles.speciesdetails_current); % De-select;
      clear WhichCheckbox;
    end  
  end; clear status;

function handles = OpenSpeciesDetailsWindow(hObject, eventdata, handles)
    
        % Yun Jiang - Debug
  if 0, fprintf('In function of OpenSpeciesDetailsWindow\n');end;
  
  
  for n=1:numel(handles.speciesdetails),
    set(handles.speciesdetails(n),'visible','off');
  end; clear n;

  if isequal(get(handles.uipanel_speciesdetails,'visible'),'off')
    % Open animation
    TotalAnimationDuration = 0.3;        % In seconds
    NumAnimationSteps = 5;
    FinalPos = [33,0.769,43,19.53846153846154];
    WhichCheckbox = handles.CurrentSpecies - handles.library_toppos + 1; % handles.CurrentSpecies is in the new entries
    if WhichCheckbox<1 | WhichCheckbox>numel(handles.checkboxes),
      WhichCheckbox = bitshift(numel(handles.checkboxes),-1)+1;
    end
    tmppos1 = get(handles.checkboxes(WhichCheckbox),'position'); clear WhichCheckbox;
    tmppos0 = get(handles.uipanel_library,'position');
    InitialPos = [tmppos0(1)+tmppos0(3), tmppos0(2)+tmppos1(2), 5, tmppos1(4)]; clear tmppos0 tmppos1;
    set(handles.uipanel_speciesdetails,'visible', 'off','Title',' ');
    tic
    AnimationStepDuration = TotalAnimationDuration/NumAnimationSteps;
    for n=1:NumAnimationSteps
      set(handles.uipanel_speciesdetails,'position',InitialPos + (FinalPos-InitialPos)/(NumAnimationSteps-1)*(n-1), 'visible', 'on');
      tmpval = n*AnimationStepDuration - toc;
      if tmpval>0, pause(tmpval); else n = ceil(toc/AnimationStepDuration); end
    end; clear n tmpval InitialPos FinalPos TotalAnimationDuration NumAnimationSteps AnimationStepDuration;
    set(handles.uipanel_speciesdetails,'Title','Details');
    for n=1:numel(handles.speciesdetails),
      set(handles.speciesdetails(n),'visible','on');
    end; clear n;
  else
    % Refresh animation
    TotalAnimationDuration = 0.3;        % In seconds
    for n=1:numel(handles.speciesdetails),
      set(handles.speciesdetails(n),'visible','off');
    end; clear n;
    pause(TotalAnimationDuration);
    for n=1:numel(handles.speciesdetails),
      set(handles.speciesdetails(n),'visible','on');
    end; clear n;
  end
  if handles.CurrentSpecies <= numel(handles.species),
    % Existing entry
    set(handles.pushbutton_deletespecies,'enable','on','visible','on');
    tmpstr = handles.species(handles.CurrentSpecies).name; set(handles.edit_speciesname,'String',tmpstr,'UserData',tmpstr);
    tmpstr = sprintf('%g\n',handles.species(handles.CurrentSpecies).frequency); set(handles.edit_ppm   ,'String',tmpstr,'UserData',tmpstr);
    tmpstr = sprintf('%g\n',handles.species(handles.CurrentSpecies).relAmps); set(handles.edit_relAmps ,'String',tmpstr,'UserData',tmpstr);
    clear tmpstr;
  else
    set(handles.pushbutton_deletespecies,'enable','off','visible','off');
    tmpstr = ''; set(handles.edit_speciesname,'String',tmpstr,'UserData',tmpstr);
    tmpstr = ''; set(handles.edit_ppm        ,'String',tmpstr,'UserData',tmpstr);
    tmpstr = ''; set(handles.edit_relAmps    ,'String',tmpstr,'UserData',tmpstr);
    clear tmpstr;
  end
  handles.speciesdetails_current = handles.CurrentSpecies;
  set(handles.pushbutton_cancelspecies,'enable','on');
  set(handles.pushbutton_savespecies  ,'enable','on');
  set(handles.text_speciesname,'Foreground',[0,0,0],'FontWeight','normal');
  set(handles.text_ppm,'Foreground',[0,0,0],'FontWeight','normal');
  set(handles.text_relAmps,'Foreground',[0,0,0],'FontWeight','normal');
  set(handles.text_status,'String','');
  handles.speciesdetails_changed = 0;

function handles = CloseSpeciesDetailsWindow(hObject, eventdata, handles)
    
            % Yun Jiang - Debug
  if 0, fprintf('In function of CloseSpeciesDetailsWindow\n');end;
  
  if ~isempty(handles.speciesdetails_current),
    for n=1:numel(handles.speciesdetails),
      set(handles.speciesdetails(n),'visible','off');
    end; clear n;

    for n=1:numel(handles.speciesdetails),
      set(handles.speciesdetails(n),'visible','off');
    end; clear n;

    % Animation
    TotalAnimationDuration = 0.3;        % In seconds
    NumAnimationSteps = 5;
    InitialPos = [33,0.769,43,19.53846153846154];
    FinalSize = [5,5];
    AnimationStepDuration = TotalAnimationDuration/NumAnimationSteps;
    FinalPos = [InitialPos(1)+(InitialPos(3)-FinalSize(1))/2,InitialPos(2)+(InitialPos(4)-FinalSize(2))/2,FinalSize(1),FinalSize(2)]; clear FinalSize;
    set(handles.uipanel_speciesdetails,'Title',' ');
    tic
    for n=1:NumAnimationSteps
      set(handles.uipanel_speciesdetails,'position',InitialPos + (FinalPos-InitialPos)/(NumAnimationSteps-1)*(n-1), 'visible', 'on');
      tmpval = n*AnimationStepDuration - toc;
      if tmpval>0, pause(tmpval); else n = ceil(toc/AnimationStepDuration); end
    end; clear n tmpval InitialPos FinalPos TotalAnimationDuration NumAnimationSteps AnimationStepDuration;
    set(handles.uipanel_speciesdetails,'visible','off');
  
    handles.speciesdetails_current = [];
    handles.speciesdetails_changed = 0;
  end

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

         % Yun Jiang - Debug
  if 0, fprintf('In function of slider1_Callback\n');end;
  
  currentpos = round(get(hObject,'Value'));
  maxval = get(hObject,'Max');
  set(hObject,'Value',currentpos);
  handles.library_toppos = maxval+1-currentpos;
  guidata(hObject, handles); % Update handles structure  
  UpdateLibraryList(hObject, eventdata, handles);
  clear currentpos maxval
  
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
         % Yun Jiang - Debug
  if 0, fprintf('In function of slider1_CreateFcn\n');end;

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_savespecies.
function pushbutton_savespecies_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_savespecies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

         % Yun Jiang - Debug
  if 0, fprintf('In function of pushbutton_savespecies_Callback\n');end;
  
  if handles.speciesdetails_changed == 0,
    set(handles.text_status,'String','Nothing to save','ForegroundColor',[0,0,0],'FontWeight','normal');
    return;
  end
  CanSave = CanSaveSpeciesDetails(hObject,eventdata,handles);
  if CanSave==0,
    return;
  end
  [status,handles] = SaveSpeciesDetails(hObject, eventdata, handles);
  clear status CanSave;
  guidata(hObject, handles);
  
% --- Executes on button press in pushbutton_cancelspecies.
function pushbutton_cancelspecies_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancelspecies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        % Yun Jiang - Debug
  if 0, fprintf('In function of pushbutton_cancelspecies_Callback\n');end;

if handles.speciesdetails_changed,
  tmpstr = strtrim(get(handles.edit_speciesname,'String'));
  if ~isempty(tmpstr), tmpstr = sprintf('in ''%s''',tmpstr); else tmpstr = 'on the right'; end
  if ~isequal(questdlg(sprintf('Discard changes %s?',tmpstr), 'Confirm',  'Cancel','Discard','Cancel'),'Discard'),
    clear tmpstr;
    return; 
  end
  clear tmpstr;
end
if handles.speciesdetails_current == numel(handles.species)+1  % New entry
  WhichCheckbox = handles.speciesdetails_current - handles.library_toppos + 1; % handles.CurrentSpecies is in the new entries
  if WhichCheckbox>=1 & WhichCheckbox<=numel(handles.checkboxes),
    set(handles.checkboxes(WhichCheckbox),'Value',0);
  end
  handles.library_selection = setdiff(handles.library_selection, handles.speciesdetails_current);   % De-select
  clear WhichCheckbox;
end
handles = CloseSpeciesDetailsWindow(hObject, eventdata, handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton_resetlibrary.
function pushbutton_resetlibrary_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resetlibrary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

       % Yun Jiang - Debug
  if 0, fprintf('In function of pushbutton_resetlibrary_Callback\n');end;

    status = questdlg('Reset library? This CANNOT be undone.', 'Confirm',  'Reset','Cancel','Cancel');
    if isequal(status,'Reset'),
      species = DefaultLibrary();
      status = WriteLibrary(handles.LibraryFilename, species);
      if status,
        errordlg({'Cannot save library file',sprintf('%s',handles.LibraryFilename)},'Error','modal');
      end
      handles = PopulateLibraryList(hObject, eventdata, handles, species);
      guidata(hObject, handles); % Update handles structure  
      clear species;
    end
    clear status;
    

