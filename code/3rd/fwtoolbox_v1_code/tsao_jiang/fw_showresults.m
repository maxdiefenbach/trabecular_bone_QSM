function fw_showresults(handles, varargin)
  if ~isfield(handles,'figure'), return; end
  if isfield(handles,'pushbutton_save') & ishandle(handles.pushbutton_save),
    set(handles.pushbutton_save,'enable','off');
  end
  if isfield(handles,'pushbutton_close') & ishandle(handles.pushbutton_close),
    set(handles.pushbutton_close,'enable','off');
  end
  if isfield(handles,'pushbutton_uponeslice') & ishandle(handles.pushbutton_uponeslice),
    set(handles.pushbutton_uponeslice,'enable','off');
  end
  if isfield(handles,'pushbutton_downoneslice') & ishandle(handles.pushbutton_downoneslice),
    set(handles.pushbutton_downoneslice,'enable','off');
  end
  if isfield(handles,'pushbutton_grayscale') & ishandle(handles.pushbutton_grayscale),
    set(handles.pushbutton_grayscale,'enable','off');
  end
  if isfield(handles,'pushbutton_colorscale') & ishandle(handles.pushbutton_colorscale),
    set(handles.pushbutton_colorscale,'enable','off');
  end

  try
    % Set figure
    if ~ishandle(handles.figure), figure(handles.figure); end
    set(handles.figure,'Pointer', 'watch');                % Wait pointer
    if isfield(handles,'DataName'),
      set(handles.figure,'name',sprintf('Click to change slice (%s)',handles.DataName));
    else
      set(handles.figure,'name','Click to change slice');
    end

    % Get information about slices
    if isfield(handles.outParams,'water'),
      NumSlices = size(handles.outParams.water,3);
    else
      NumSlices = size(handles.outParams.species(1).amps,3);
    end
    if ~isfield(handles,'z'),
      handles.z = bitshift(NumSlices,-1)+1;
    end

    % Buttons
    x=0;
    if ~isfield(handles,'pushbutton_downoneslice') | ~ishandle(handles.pushbutton_downoneslice),
      width = 5;
      handles.pushbutton_downoneslice = uicontrol(handles.figure,'Style','Pushbutton','Units','Characters','position',[x,0,width,1.5],'String','<', ...
          'Callback','fw_showresults(guidata(gcbo),''downoneslice'')','enable','off'); 
      x=x+width;
    end
    if ~isfield(handles,'pushbutton_uponeslice') | ~ishandle(handles.pushbutton_uponeslice),
      width = 5;
      handles.pushbutton_uponeslice = uicontrol(handles.figure,'Style','Pushbutton','Units','Characters','position',[x,0,width,1.5],'String','>',...
          'Callback','fw_showresults(guidata(gcbo),''uponeslice'')','enable','off'); 
      x=x+width;
    end
    if ~isfield(handles,'pushbutton_grayscale') | ~ishandle(handles.pushbutton_grayscale),
      width = 10;
      handles.pushbutton_grayscale = uicontrol(handles.figure,'Style','Pushbutton','Units','Characters','position',[x,0,width,1.5],'String','Gray',...
          'Callback','fw_showresults(guidata(gcbo),''grayscale'')','enable','off'); 
      x=x+width;
    end
    if ~isfield(handles,'pushbutton_colorscale') | ~ishandle(handles.pushbutton_colorscale),
      width = 10;
      handles.pushbutton_colorscale = uicontrol(handles.figure,'Style','Pushbutton','Units','Characters','position',[x,0,width,1.5],'String','Color',...
          'Callback','fw_showresults(guidata(gcbo),''colorscale'')','enable','off'); 
      x=x+width;
    end
    if ~isfield(handles,'pushbutton_save') | ~ishandle(handles.pushbutton_save),
      width = 10;
      handles.pushbutton_save = uicontrol(handles.figure,'Style','Pushbutton','Units','Characters','position',[x,0,width,1.5],'String','Save...', ...
          'Callback','fw_showresults(guidata(gcbo),''save'')','enable','off'); 
      x=x+width;
    end
    if ~isfield(handles,'pushbutton_close') | ~ishandle(handles.pushbutton_close),
      width = 10;
      handles.pushbutton_close = uicontrol(handles.figure,'Style','Pushbutton','Units','Characters','position',[x,0,width,1.5],'String','Close',...
          'Callback','fw_showresults(guidata(gcbo),''close'')','enable','off'); 
      x=x+width;
    end
    clear x width;    
    
    % Action
    if ~isempty(varargin),
      if ischar(varargin{1}),
        switch lower(varargin{1}),
          case 'mouse'
            % Check changes in slice number
            WhichMouseButton = get(handles.figure, 'SelectionType');
            if isequal(WhichMouseButton,'alt'),                       % Right-clicked in images
              handles.z=handles.z-1;
              if handles.z<1, handles.z=NumSlices; end
            elseif isequal(WhichMouseButton,'normal'),                % Left-clicked in images
              handles.z=handles.z+1;
              if handles.z>NumSlices, handles.z=1; end
            end
            clear WhichMouseButton;
          case 'uponeslice'
            handles.z=handles.z+1;
            if handles.z>NumSlices, handles.z=1; end
          case 'downoneslice'
            handles.z=handles.z-1;
            if handles.z<1, handles.z=NumSlices; end
          case 'grayscale'
            colormap(gray);
          case 'colorscale'
            colormap(jet);
          case 'close'
            delete(handles.figure); return;
          case 'save'
            str1 = 'Save to console';
            str2 = 'Save to .mat file';
            ButtonName=questdlg('Where to save?', ...
                       'Saving...', str1, str2, str1);
            if isequal(ButtonName,str1)
              assignin('base','outParams',handles.outParams);
              if isfield(handles,'outParamsMP'),
                assignin('base','outParamsMP',handles.outParamsMP);
              end
            elseif isequal(ButtonName,str2)
              [tmpfile, tmppath] = uiputfile({'*.mat', 'MATLAB File (*.mat)'}, 'Save as');
              if ~isequal(tmpfile,0)
                if isfield(handles,'outParamsMP'),
                  save(fullfile(tmppath,tmpfile), '-STRUCT', 'handles', 'outParams', 'outParamsMP');
                else
                  save(fullfile(tmppath,tmpfile), '-STRUCT', 'handles', 'outParams');
                end
              end
              clear tmpfile tmppath;
            end
            clear str1 str2;
        end
      end
    end

    % Set up figure
    if isfield(handles.outParams,'water'),
      % Original algorithm with single-peak and multi-peak results
      h=subplot(2,4,1,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.water(:,:,handles.z))','Parent',h);title(h,{'WATER image',sprintf('z = %d',handles.z)},'vertical','bottom'); ylabel(h,'Single-peak'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      h=subplot(2,4,2,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.fat(:,:,handles.z))','Parent',h);title(h,{'FAT image',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      h=subplot(2,4,3,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.r2starmap(:,:,handles.z))','Parent',h);title(h,{'R2*',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      h=subplot(2,4,4,'Replace',handles.figure); h2=imagesc(angle(handles.outParams.phasemap(:,:,handles.z))','Parent',h);title(h,{'Phase map',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      if isfield(handles,'outParamsMP'),
        h=subplot(2,4,5,'Replace',handles.figure); h2=imagesc(abs(handles.outParamsMP.species(1).amps(:,:,handles.z)'));title(h,{'WATER image',sprintf('z = %d',handles.z)},'vertical','bottom'); ylabel(h,'Multi-peak'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
        set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
        h=subplot(2,4,6,'Replace',handles.figure); h2=imagesc(abs(handles.outParamsMP.species(2).amps(:,:,handles.z)'));title(h,{'FAT image',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
        set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
        h=subplot(2,4,7,'Replace',handles.figure);
        h2=imagesc(abs(handles.outParamsMP.SPtoMPmatrix),'Parent',h); axis(h,'image');
        set(h,'XTick',[1,2,3],'XTickLabel',{'Water','Fat','Blank'},'YTick',[1,2],'YTickLabel',{'Water','Fat'});
        xlabel(h,'Single peak');
        ylabel(h,'Multi-peak');
        title(h,'Conversion matrix');
        set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      end
      drawnow;
      h=subplot(2,4,8,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.fiterror(:,:,handles.z))','Parent',h);title(h,{'Fitting Error',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
    else
      % New algorithm with multi-peak results only
      for n=1:numel(handles.outParams.species),
        h=subplot(2,numel(handles.outParams.species),n,'Replace',handles.figure); 
        h2=imagesc(abs(handles.outParams.species(n).amps(:,:,handles.z)'));
        title(h,{handles.outParams.species(n).name,sprintf('z = %d',handles.z)},'vertical','bottom'); 
        ylabel(h,''); axis(h,'image'); 
        h3=colorbar('peer',h); drawnow;
        set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); 
        set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      end; clear n;
%      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
%      h=subplot(2,2,2,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.species(2).amps(:,:,handles.z)'));title(h,{'FAT image',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
%      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      h=subplot(2,3,4,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.r2starmap(:,:,handles.z))','Parent',h);title(h,{'R2*',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      h=subplot(2,3,5,'Replace',handles.figure); h2=imagesc(angle(handles.outParams.phasemap(:,:,handles.z))','Parent',h);title(h,{'Phase map',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
      h=subplot(2,3,6,'Replace',handles.figure); h2=imagesc(abs(handles.outParams.fiterror(:,:,handles.z))','Parent',h);title(h,{'Fitting Error',sprintf('z = %d',handles.z)},'vertical','bottom'); h3=colorbar('peer',h); axis(h,'image'); drawnow;
      set(h2,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')'); set(h3,'ButtonDownFcn','fw_showresults(guidata(gcbo),''mouse'')');
    end
    clear h h2 h3;
  end    
  
  if isfield(handles,'pushbutton_save') & ishandle(handles.pushbutton_save),
    set(handles.pushbutton_save,'enable','on');
  end
  if isfield(handles,'pushbutton_close') & ishandle(handles.pushbutton_close),
    set(handles.pushbutton_close,'enable','on');
  end
  if isfield(handles,'pushbutton_uponeslice') & ishandle(handles.pushbutton_uponeslice),
    set(handles.pushbutton_uponeslice,'enable','on');
  end
  if isfield(handles,'pushbutton_downoneslice') & ishandle(handles.pushbutton_downoneslice),
    set(handles.pushbutton_downoneslice,'enable','on');
  end  
  if isfield(handles,'pushbutton_grayscale') & ishandle(handles.pushbutton_grayscale),
    set(handles.pushbutton_grayscale,'enable','on');
  end
  if isfield(handles,'pushbutton_colorscale') & ishandle(handles.pushbutton_colorscale),
    set(handles.pushbutton_colorscale,'enable','on');
  end
  if ishandle(handles.figure), 
    set(handles.figure,'Pointer', 'arrow', 'visible','on');                % Arrow pointer 
    guidata(handles.figure, handles);
  end  
  drawnow;  
