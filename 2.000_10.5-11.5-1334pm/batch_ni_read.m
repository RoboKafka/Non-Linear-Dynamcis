 function varargout=batch_ni_read(varargin)
%
% NI signal reader
%

% BT/UOB/22102013, Last revision: 29/01/2014

happ=findobj('tag','batch_ni_read_figure');

if nargin==0
   % initialization
   if isempty(happ)
      action='init';
   else
      error('batch_ni_read is already open.')
   end
elseif nargin==1
   % direct calls
   if strcmp(varargin{1},'initlist') || strcmp(varargin{1},'init')
      action=varargin{1};
   else
      error('Incorrect direct call.')
   end
else
   % callback-styled calls
   action=varargin{3};
   if ~isempty(happ)
      ud=get(happ,'userdata');
   else
      error('NI read application not available.')
   end
end

switch action
   case 'exec'
      % let's assume this function exists within, no try-catch!
      datout=cell(1,nargout);
      [datout{:}]=feval(varargin{4:end});
      varargout=datout;return
        
   case 'init'
      if ~isempty(happ)
         error('batch_ni_read is already open.')
      end
      l_init
      l_window
      
   case 'initlist'
      batch_ni_read('init')
      happ=findobj('tag','batch_ni_read_figure');
      l_listbox('batch_ni_read_figure_list',happ);
      
   case {'start','stop'}
      ses=ud.ni_session_obj;
      if isempty(ses)
         error('NI session not established.')
      end
      % switching action bold on buttons
      switch action
         case 'start'
            if ~ud.is_run
               ses.startBackground();
               ud.is_run=true;
            else
               fprintf('\nNI acquisition already running.\n')
            end
         case 'stop'
            if ud.is_run
               ses.stop();
               ud.is_run=false;
            else
               fprintf('\nNI acquisition already stopped.\n')
            end
      end
      % update app structure
      set(ud.handles{1},'userdata',ud)
      
   case 'initni'
      % check presence of NI session
      ses=ud.ni_session_obj;
      if ~isempty(ses) disp('NI session already established.'),return,end
      
      % UI section
      prompt={'Define channels (0,1,2,3):'};
      name='Input session parameters'; numlines=1;
      defaultanswer={'0'};
      answer=inputdlg(prompt,name,numlines,defaultanswer);
      if isempty(answer) return, end
      if isempty(answer{1})
         warning('No channel info provided. No NI session established.')
         return
      else
         % no checking (valid integers only!)
         chans=str2num(answer{1});
      end
      
      % -----------------------------
      %     NI instrument calls
      % -----------------------------
      % ... find/access NI instrument
      devs=daq.getDevices;       % read all devices
      dev=devs(1);               % assumes one instrument only
      % ... ... error trap here
      disp('NI read: NI instrument found')
      % ... open NI session
      ses=daq.createSession('ni');
      % ... ... error trap here
      disp('NI read: NI session established')
      % ... open channels
      for ichan=1:length(chans)
         ses.addAnalogInputChannel(dev.ID,chans(ichan),'IEPE');
         ud.axb_cal{ichan}=[1,0];
         % ... ... error trap here
         fprintf('\nNI read: ... channel %d open',chans(ichan))
      end
      ud.chan_id=chans;
      % ... generic NI session auto callback event
      listenses=ses.addlistener('DataAvailable',@l_event);
      % ... ... error trap here
      
      % NI session properties
      ses=l_set_props(ses, ...
         'IsContinuous',true, ...
         'NotifyWhenDataAvailableExceeds',ud.sa_evnt, ...
         'Rate',ud.sa_rate);
      % store NI objects
      ud.ni_device_obj=dev;
      ud.ni_session_obj=ses;
      ud.ni_listener_obj=listenses;
      % update application data
      set(ud.handles{1},'userdata',ud)
      % identify session status & deactivate callback
      fcn=get(ud.handles{3}(3),'buttondownfcn');
      set(ud.handles{3}(3), ...
         'fontweight','bold','color',[1,1,1]*0.45, ...
         'userdata',fcn,'buttondownfcn',[])
      
   case 'props'
      % get session info
      ses=ud.ni_session_obj;
      if isempty(ses) error('NI session not established.'),end
      if ud.is_run fprintf('\nChange NI properties only when stopped.'), return,end
      chans=ud.chan_id;
      
      % UI section
      prompt={ ...
         'Define *minimum* time frame in [s]:', ...
         'Define read rate (sampling) in [Sa/s]:', ...
         'Define refresh rate [Sa/s]:'};
      name='Input session parameters'; numlines=1;
      defaultanswer={ ...
         num2str(ud.dt_read) , ...
         num2str(ud.sa_rate) , ...
         num2str(ud.sa_evnt)};
      % augment nchan calib etc. lines
      for ii=1:length(chans)
         prompt=[prompt,{sprintf( ...
            'Channel %d: (x+A)*B calib ([EU/V],[EU]):',chans(ii))}];
         defaultanswer=[defaultanswer,{num2str(ud.axb_cal{ii})}];
      end
      % choose callback functions
      % ... line click callback
      prompt=[prompt,{'Line click callback (1-fft,2-frf):'}];
      defaultanswer=[defaultanswer,{num2str(ud.fcn_lineclickfcnid)}];
      % ... trigger callback @ properties [0|1,threshold]
      prompt=[prompt,{'Trigger spec ([flag_01,channel,level]):'}];
      defaultanswer=[defaultanswer,{sprintf('%d, %d, %f', ...
          ud.trig_flag,ud.trig_chan,ud.trig_lev)}];
      % ... event callback function
      prompt=[prompt,{'Event callback function (0-none,1-swft):'}];
      defaultanswer=[defaultanswer,{num2str(ud.fcn_evntcallfcnid)}];
      
      % UI PANNEL
      % ---------
      answer=inputdlg(prompt,name,numlines,defaultanswer);
      if isempty(answer) return, end
      % no further logical checks
      
      % Read UI PANNEL
      % ... time frame
      if isempty(answer{1}) warning('No time frame change.')
      else ud.dt_read=str2num(answer{1}); end
      % ... sampling rate
      if isempty(answer{2}) warning('No read rate change.')
      else
         ud.sa_rate=str2num(answer{2});
         ses=l_set_props(ses,'Rate',ud.sa_rate);
         %set(ses,'Rate',ud.sa_rate)
      end
      % ... event rate in # Sa
      if isempty(answer{3}) warning('No refresh rate change.')
      else
         ud.sa_evnt=str2num(answer{3});
         ses=l_set_props(ses,'NotifyWhenDataAvailableExceeds',ud.sa_evnt);
         %set(ses,'NotifyWhenDataAvailableExceeds',ud.sa_evnt)
      end
      % ... calibration constants
      for ii=1:length(chans)
         if isempty(answer{3+ii})
            warning('No calibration change.')
         else
            ud.axb_cal{ii}=str2num(answer{3+ii});
            fprintf('\nNI Channel %d AxB is [%.2f,%.2f]', ...
               chans(ii),ud.axb_cal{ii}(1),ud.axb_cal{ii}(2))
         end
      end
      nansw=3+length(chans);
      % ... line click callback
      if isempty(answer{1}) warning('No line click callback fcn ID.')
      else ud.fcn_lineclickfcnid=round(str2num(answer{nansw+1})); end
      fprintf('\nLine click trigger fcn=%d',ud.fcn_lineclickfcnid)
      % ... trigger spec
      if isempty(answer{1}) warning('No trigger spec provided.')
      else
          trigspecs=str2num(answer{nansw+2});
          ud.trig_flag=logical(trigspecs(1));
          ud.trig_chan=round(trigspecs(2));
          ud.trig_lev=trigspecs(3);
          fprintf('\nTrigger is *%d* on chan %d, treshold=%f', ...
             ud.trig_flag,ud.trig_chan,ud.trig_lev)
      end
      % ... event callabck spec
      if isempty(answer{1}) warning('No event callback fcn ID.')
      else ud.fcn_evntcallfcnid=round(str2num(answer{nansw+3})); end
      fprintf('\nEvent callback fcn=%d',ud.fcn_evntcallfcnid)
      fprintf('\n')

      
      % **mimum read time*** frame to int multiple of refresh rate
      ud.nblocks=ceil((ud.sa_rate*ud.dt_read)/ud.sa_evnt);
      ud.dt_read=(ud.sa_evnt/ud.sa_rate)*ud.nblocks;
      
      % update user data
      ud.ni_session_obj=ses;
      set(ud.handles{1},'userdata',ud)
      % update visuals
      l_window
      
      % Qs:
      % + can we change session pars on flight
      % + can we pass object between aliases
      
   case 'export'
      l_export_data
      
   case 'closeni'
      % close NI session
      if isempty(ud.ni_session_obj)
         fprintf('\n NI session already closed!\n')
         return
      end
      fprintf('\n\n Closing NI session!')
      delete(ud.ni_session_obj)
      % ... error trap here
      
      % clean up and reset
      ud=l_init_ud(ud.handles{1:3},ud.handles{5});
      set(ud.handles{1},'userdata',ud)
      l_window
      fcn=get(ud.handles{3}(3),'userdata');
      set(ud.handles{3}(3),'fontweight','normal','color','k', ...
         'userdata',[],'buttondownfcn',fcn)
      fprintf('\n Done!\n')
      
   otherwise
      error('Unknown call')
end
return


% --------------------------------------------
% Todo: unitwin, expwin, H1, H2, etc. FRF est, COH
%
% Main:
%  + exec
%  + init
%  + initlist
%  + start
%  + stop
%  + initni
%  + props
%  + export
%  + closeni
%  + (unknown)
%
% Subfunctions:
%  l_export_data
%
%  l_proc_lst
%  l_proc_eswft
%  l_proc_trig
%  l_proc_fft
%  l_proc_frf
%
%  l_comp_freqan    MAIN DSP FCN
%
%  l_set_props
%  l_event          MAIN NI IO FCN
%
%  l_listbox        APP4 LIST UI
%  l_window         APP3 RESTART
%  l_init           APP1 UI
%  l_init_ud        APP2 DATA
% --------------------------------------------

function l_export_data
% Export available data to workspace

happ=findobj('tag','batch_ni_read_figure');
ud=get(happ,'userdata');

dataout=[ud.time,ud.data];
dataout(isnan(ud.time),:)=[];

assignin('base','nirawdata',dataout)
fprintf('\nNI raw signals exported to workspace in *nirawdata*.\n')
return
% --------------------------------------------

function l_proc_lst(in1,in2,in3)
% Process stored data (list object db)

hl=get(get(in1,'parent'),'userdata');
val=get(hl,'value');
str=get(hl,'string');
ud=get(hl,'userdata');

% consistency pre-checks
switch in3
   case {'showsel','showselavg','expselavg','addnote','delsel'}
      if isempty(val),warning('No list item selected!'),return,end
      if isempty(str),warning('No list data available!'),return,end
end

% listbox contextmenu functionality
switch in3
   case {'showsel','showselavg'}
      % check for extra functionality
      IS_FRFUI=false;if exist('plot_uifrf','file')==2,IS_FRFUI=true;end

      % show all frfs for blocks
      figure, ha1=axes; hold on, box on, grid on
      
      for ii=1:length(val)
         db=ud.dblock{val(ii)}; col=rand(1,3); Hii=[];
         for jj=1:length(db.trig_x)
            %plot(db.trig_x{jj},db.trig_y{jj},'color',col)
            [freq,hacc]=l_comp_freqan('rawfrf', ...
               db.trig_x{jj}.',db.trig_y{jj}(:,[2,1]).');
            % output
            fcnstr=sprintf('disp(''%s @ %s'')', ...
               str{val(ii)},ud.timestamp{val(ii)});
            switch in3
               case 'showsel'
                  hl=plot(freq(1:floor(end/2)), ...
                     log10(abs(hacc(1:floor(end/2)))),'color',col, ...
                     'buttondownfcn',fcnstr);
                  % add FRF ui capability
                  if IS_FRFUI,set(hl,'userdata',{freq,hacc}),end
               case 'showselavg'
                  Hii=[Hii;hacc];
            end
         end
         % show mean values if requested
         if strcmp(in3,'showselavg')
            Hiim=mean(Hii,1);
            hl=plot(freq(1:floor(end/2)), ...
               log10(abs(Hiim(1:floor(end/2)))),'color',col, ...
               'buttondownfcn',fcnstr);
            % add FRF ui capability
            if IS_FRFUI,set(hl,'userdata',{freq,hacc}),end
         end
      end
      % add FRF ui capability
      if IS_FRFUI, plot_uifrf(ha1,[],'initui'), end

   case 'cfitone'
      % circle fitting in external function
      for ii=1:length(val)
         db=ud.dblock{val(ii)}; col=rand(1,3); Hii=[];
         for jj=1:length(db.trig_x)
            % recalculate FRF
            [freq,hacc]=l_comp_freqan('rawfrf', ...
               db.trig_x{jj}.',db.trig_y{jj}(:,[2,1]).');
            % output
            fcnstr=sprintf('disp(''%s @ %s'')', ...
               str{val(ii)},ud.timestamp{val(ii)});
            % CFIT routine link
            freq1=freq(:);hvel=hacc(:)./(sqrt(-1)*2*pi*freq1);
            plot_cfitone([freq1,hvel])
         end
      end
      
   case 'expselavg'
      disp('export selected mean frf etc to Matlab ws')
      
   case 'addnote'
      % additional data block description
      note=inputdlg('Add a note to all selected list items:','Note',1,{''});
      if isempty(note) return, end
      if ~isempty(note{1})
         for ii=1:length(val)
            str{val(ii)}=sprintf('%d: %s',val(ii),note{1});
            ud.textnote{val(ii)}=note{1};
         end
         set(hl,'string',str)
      end
      set(hl,'userdata',ud); % isn't this op too expensive?
      
   case 'saveud'
      % save acquired listed data
      %uisave('ud',sprintf('nidata_%s.nid',datestr(now,30)))
      ftitle=sprintf('nidata_%s.nid',datestr(now,30));
      [fname,pname]=uiputfile({'*.nid','NI-data files (*.nid)'}, ...
         'Save as',ftitle);
      if isequal(fname,0) || isequal(pname,0)
         return
      else
         ffname=fullfile(pname,fname); %disp(['Saved file: ',ffname])
         save(ffname,'ud','-mat')
      end
      
   case 'loadud'
      % read saved list data
      [filename,pathname]=uigetfile('*.nid', 'Pick a NI-data file');
      if isequal(filename,0) || isequal(pathname,0)
         return
      end
      
      try
         disp(['Reading: ',fullfile(pathname,filename)])
         data=load(filename,'-mat');
         %fnm1=fieldnames(data);fnm2=fieldnames(data.(fnm1{1}));
         udnew=data.('ud');    % assumes this is expected ud!
         
         for ii=1:length(udnew.dbid)
            str{ii}=sprintf('%d: %s',udnew.dbid(ii),udnew.textnote{ii});
         end
         set(hl,'userdata',udnew,'string',str)
         
      catch
         disp('ooops error, can not process selected file')
         return
      end
      
      % retain this for new 'conditional' uicontext menu item!!!
      if false
      for ii=1:length(val)
         db=ud.dblock{val(ii)}; col=rand(1,3); Hii=[];
         for jj=1:length(db.trig_x)
            %plot(db.trig_x{jj},db.trig_y{jj},'color',col)
            read_freevibproc( ...
               [db.trig_x{jj}(:),db.trig_y{jj}(:,[1])]);
            % output
            %fcnstr=sprintf('disp(''%s @ %s'')',str{val(ii)},ud.timestamp{val(ii)});
         end
      end
      end
      
   case 'delsel'
      % delete items
      dlt=questdlg('Delete selected items?','Delete','Ok','Cancel','Cancel');
      if isempty(dlt) || strcmp(dlt,'Cancel'), return, end
      % apply delete (counter is universal)
      ud.dbid(val)=[];
      ud.timestamp(val)=[];
      ud.dblock(val)=[];
      str(val)=[];
      set(hl,'string',str,'userdata',ud,'value',[])
end

return
% --------------------------------------------

function l_proc_eswft(xx,yy)
% Event callback function: SWFT
% not UI fcn, triggered with event

figname='batch_ni_read_figure_swft';
hfig=findobj('tag',figname);
if isempty(hfig)
   hfig=figure('tag',figname);ha=axes('ylim',[-4,1]);
   hold on,box on,grid on
   xlabel('frequency [hz]'), ylabel('log10(abs(chan))')
   hl=line('xdata',[],'ydata',[],'parent',ha);
   set(hfig,'userdata',hl)
else
   hl=get(hfig,'userdata');
end

figname2='batch_ni_read_figure_swft2';
hfig2=findobj('tag',figname2);
if isempty(hfig2)
   hfig2=figure('tag',figname2);ha2=axes;
   hold on,box on,grid on
   set(hfig2,'userdata',ha2)
else
   ha2=get(hfig2,'userdata');
end

% do something here
[ff,YY]=l_comp_freqan('fft',xx',yy(:,1)'); ff=ff(:); YY=YY.';
% basic plotting (new figure)
set(hl, ...
   'xdata',ff(1:floor(end/2)), ...
   'ydata',log10(abs(YY(1:floor(end/2)))))
drawnow

% new style test
nfp=floor(length(ff)/2);
hsf=surf(ha2, ...
   ones(nfp,1) * [xx(1),xx(1)+(xx(50)-xx(1))], ...
   ff(1:nfp)             * [1,1], ...
   log10(abs(YY(1:nfp))) * [1,1], ...
   log10(abs(YY(1:nfp))) );
set(hsf,'edgecolor','none')
drawnow
return
% --------------------------------------------

function l_proc_trig(xtrig,ytrig,flag)
% process trigger data - TEST IMPLEMENTATIONS

figname='batch_ni_read_figure_trigdat';
appname='batch_ni_read_figure';
lisname='batch_ni_read_figure_list';
hfig=findobj('tag',figname);
hfig2=findobj('tag',appname);
hobj=findobj('tag',lisname);
udapp=get(hfig2,'userdata');

if isempty(hfig)
   % create figure for one acquisition
   hfig=figure('tag',figname,'WindowKeyPressFcn',{@l_proc_trig,'key'}); %KeyPressFcn
   ha=axes;
   xy0=get(findobj('tag','batch_ni_read_figure'),'position');
   xy1=get(hfig,'position');
   set(hfig,'position',[xy0(1)+xy0(3),xy0(2),xy1(3:4)])
   ud=struct('trig_x',{{}},'trig_y',{{}});  % _trigdat ud template
   set(hfig,'userdata',ud)
else
   % process call after trigger window
   ud=get(hfig,'userdata');
   ha=get(hfig,'children');
   % process key press
   if strcmp(flag,'key')
      current_key=get(hfig,'CurrentCharacter');
      if current_key=='n'
         % reject last stored record (if any)
         n_dblocks=length(ud.trig_x);
         if n_dblocks==0
            disp('No data in averager ...')   % should not happen
         elseif n_dblocks==1
            delete(hfig)     % kill the figure if no more data to delete
         else
            % remove last & emulate run after acquisition
            xtrig=ud.trig_x{end-1}; ud.trig_x([end-1,end])=[];
            ytrig=ud.trig_y{end-1}; ud.trig_y([end-1,end])=[];
            set(hfig,'userdata',ud)
            l_proc_trig(xtrig,ytrig,'run')
         end
      elseif current_key=='s'
         % store acquired data
         if isempty(hobj)
            [udlis,hobj]=l_listbox(lisname,hfig2); % create list object
         else
            udlis=get(hobj,'userdata'); % get listobject data for update
         end
         % update data block container
         udlis.counter=udlis.counter+1;
         udlis.dbid(end+1)=udlis.counter;
         udlis.timestamp{end+1}=datestr(now,30);
         udlis.textnote{end+1}='';
         udlis.dblock{end+1}=ud;
         % update list object
         set(hobj,'string',[get(hobj,'string');{num2str(udlis.counter)}])
         set(hobj,'userdata',udlis)
         
      end
      return % end this exec branch
   end
end

% augment new data
ud.trig_x{end+1}=xtrig;
ud.trig_y{end+1}=ytrig;

n_dblocks=length(ud.trig_x);
disp(sprintf('Stored %d trigger data blocks.',n_dblocks))

% reorder signals
ii_ord=[1:length(udapp.chan_id)];
ii_ham=find(udapp.chan_id==udapp.trig_chan);
ii_ord(ii_ham)=[];
ii_ord=[ii_ham,ii_ord];

% plots
delete(get(ha,'children')), hold on, HACCall=[];
for ii=1:n_dblocks
   % computer raw frf
   [ff,HACC]=l_comp_freqan('rawfrf',ud.trig_x{ii}',ud.trig_y{ii}(:,ii_ord)');
   HACCall=[HACCall;HACC];
   % Time plot
   %plot(ud.trig_x{ii}-ud.trig_x{ii}(1),ud.trig_y{ii},'color','b','parent',ha)
   % Frequency plot
   plot(ff(1:floor(end/2)),log10(abs(HACC(1:floor(end/2)))), ...
      'color',[1,1,1]*(3+rand(1)*4)/10,'parent',ha)
end
% brute force complex data stats 
HACCmean=mean(HACCall,1);              % mean value
plot(ff(1:floor(end/2)),log10(abs(HACCmean(1:floor(end/2)))), ...
   'color','b','parent',ha)

% update data structure
set(hfig,'userdata',ud)
return
% --------------------------------------------

function l_proc_fft(hcaller,in2)
% simple fft on non-NaN caller data

xx=get(hcaller,'xdata'); xx(isnan(xx))=[];
yy=get(hcaller,'ydata'); yy(isnan(yy))=[];
if isempty(xx), disp('No data to process.'),return,end

% data processing
[ff,YY]=l_comp_freqan('fftrel',xx,yy);
% data plotting (new figure)
figure, grid on, box on
plot(ff(1:floor(end/2)),log10(abs(YY(1:floor(end/2)))))
xlabel('frequency [hz]'), ylabel('log10(abs(chan))')
return
% --------------------------------------------

function l_proc_frf(hcaller,in2)
% 'FRF' on read signals, caller is hammer

happ=findobj('tag','batch_ni_read_figure');
ud=get(happ,'userdata');

% find hammer and other signals
hsigs=ud.handles{4};
iham=find(hsigs==hcaller);   % hammer channel is caller line!
if isempty(iham)
   error('error in frf callback')
else
   hham=hsigs(iham);hsigs(iham)=[];hacc=hsigs;nacc=length(hacc);
   hall=[hham,hacc];
end

% prepare data for processing
xx=get(hham,'xdata'); jnk=isnan(xx); xx(jnk)=[];
yy=[];for ii=1:length(hall),yy=[yy;get(hall(ii),'ydata')];end,yy(:,jnk)=[];

% data processing
[ff,HACC]=l_comp_freqan('rawfrf',xx,yy);
% data plotting (new figure)
figure, grid on, box on
plot(ff(1:floor(end/2)),log10(abs(HACC(:,1:floor(end/2)))))
xlabel('frequency [hz]'), ylabel('log10(abs(Y/X))')
return
% --------------------------------------------

function [xf,yh]=l_comp_freqan(faflag,tdata,sdata)
% Frequency analysis

% frequency domain discretisation
fsa=1/(tdata(2)-tdata(1));
xf=linspace(0,1,length(tdata))*fsa;

switch faflag
   case 'fft'
      yh=fft(sdata);
      
   case 'fftrel'
      frq_ref=30;  % low-pass reference freq
      yh=fft(sdata);
      yh=yh/mean(abs(yh(xf<=frq_ref)));

   case 'rawfrf'
      % FORMAT: data in rows, first row hammer
      nacc=size(sdata,1)-1;
      % first channel is hammer
      YHAM=fft(sdata(1,:));
      YACC=fft(sdata(2:end,:));
      % ... raw analytical 'frf', no psd and stats
      yh=YACC./(ones(nacc,1)*YHAM);
      
   case 'frf'
      % H1 or H2 estimator
      
end

return
% --------------------------------------------

function obj=l_set_props(obj,varargin)
% set NI session parameters (assumes correct format!)

for ii=1:(length(varargin)/2)
   % get andset prop/value pair
   prop=varargin{2*ii-1};
   value=varargin{2*ii};
   set(obj,prop,value)
   % additional info
   switch prop
      case 'IsContinuous'
         fprintf('\nNI read: ... ... continuous is: %d',true) % FIX this
      case 'NotifyWhenDataAvailableExceeds'
         fprintf('\nNI read: ... ... refresh %d Sa',value)
      case 'Rate'
         fprintf('\nNI read: ... ... read rate %d Sa/s',value)
   end
end
fprintf('\n')
return
% --------------------------------------------

function l_event(src,event)
% NI event managed by the listener

%persistent TRIGDB
%if isempty()

happ=findobj('tag','batch_ni_read_figure');
ud=get(happ,'userdata');
% ... data buffers
timebuffer=ud.time;
databuffer=ud.data;
bdatbuffer=ud.bdat;
% ... trigger parameters
IS_TRIG=ud.trig_flag;
idx_hammer_chan=ud.trig_chan+1;  % trigger channel, start counting from 0
trig_level=ud.trig_lev;          % trigger level threshold
% ... event callback
idxeventcall=ud.fcn_evntcallfcnid+1;
fcn_eventcall=ud.fcn_evntcallhans{idxeventcall}; % no call if empty

% READ EVENT OBJECT
% -----------------
xx=event.TimeStamps;
yy=event.Data;
% ... time vector prep
xx=xx-xx(1);                             % local time: 0 to tblock=sa_evnt*sa_rate
ud.t_abs=ud.t_abs+xx(end)+1/ud.sa_rate;  % global time
xxabs=xx+ud.t_abs;                       % absolute time vector

% EVENT CALLBACK
% --------------
if ~isempty(fcn_eventcall)
   fcn_eventcall(xxabs,yy) % ... no signal ordering applied!
end

% BASIC DATA HANDLING
% -------------------
ud.bcount=ud.bcount+1;
idxstart=ud.sa_end+1;
idxend=ud.sa_end+ud.sa_evnt;
if idxend<=length(timebuffer)
   ud.sa_end=idxend;
   % time vector adjustment
   xx=xx+ud.t_end;               % shift to transform to global time
else                             % reset counters to 0 in figure
   ud.bcount=1;
   ud.sa_end=0;
   idxstart=1;
   idxend=ud.sa_end+ud.sa_evnt;
   ud.sa_end=idxend;
   % time vector adjustment
   % ... global time=local time
   % clean the full buffer (no FIFO etc.)
   timebuffer(:,:)=NaN;
   databuffer(:,:)=NaN;
end
% update data content
timebuffer(idxstart:idxend,1)=xx;
databuffer(idxstart:idxend,:)=yy;
bdatbuffer{ud.bcount}=[xxabs,yy];
% future time reference
ud.t_end=xx(end)+1/ud.sa_rate;  % new refernce old + 1 sample


% SW TRIGGER HANDLING
% -------------------
% IS TRIG: TRUE, in not ON is DETECTED
if IS_TRIG && ...
      any(yy(:,idx_hammer_chan)>trig_level) && ...
      ud.trig_cnt==0
   % how to delay information about trigger to future events
   % i.e. to evaluate the trigger based on what will come next
   % e.g. multiple hits int the next block
   disp('trigger ON ...')
   ud.trig_cnt=ud.trig_cnt+1;
   ud.trig_idx=ud.bcount;
   
   % start trig counter to completion
   % and execute: stop, store, proc
   
   % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIG PLOTTING
   % plot trigger block
   pos=get(ud.handles{1},'position');
   figure('position',[pos(1),pos(2)/3,pos(3),pos(2)/2])
   hold on, box on, grid on
   plot(xxabs,yy(:,idx_hammer_chan),'b.-','linewidth',1)
   % plot pre-trigger block
   if ud.bcount>1
      xy=bdatbuffer{ud.bcount-1};
      plot(xy(:,1),xy(:,3),'m.-')
   end
   % plot post-trigger block
   if ud.bcount<ud.nblocks
      
   end
   % plot treshold line
   xlim=get(gca,'xlim');
   plot(xlim,[1,1]*trig_level,'r--')
   % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
end
% IS TRIG: TRUE, ON, READING and STOP
if IS_TRIG && ...
      ud.trig_cnt>0
   % complete trigger event
   ud.trig_cnt=ud.trig_cnt+1;
   trigtext=sprintf('%.0f%%',floor(100*ud.trig_cnt/ud.nblocks));
   set(ud.handles{5},'string',trigtext)
   % detect end of reading
   if ud.trig_cnt>=ud.nblocks%-1 % do not kill pre-trigger info
      % plot one full acquisition
      if ud.trig_idx==1
         Dc=bdatbuffer(1:end);
      else
         % ringbuffer reorder
         Dc=bdatbuffer([ud.trig_idx-1:end,1:ud.trig_idx-2]);
      end
      
      % adjust time vector here
      % ... should work with absolute time here!!!
      %
      % THIS IS CHAN DEPENDENT CORRECT THIS
      % 
      D=[Dc{:}]; nsep=3;
      % FIX THIS
      X =D(:,1:nsep:end);  % time
      Y1=D(:,2:nsep:end);  % a1
      Y2=D(:,3:nsep:end);  % f0
      %Y3=D(:,4:nsep:end);  % a2
      Y3=[];
      % FIX THIS
      XRUN=X(:);
      YRUN=[Y1(:),Y2(:),Y3(:)];
      % Note: the use of global time should ensure correct data order
       
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIG PLOTTING
      % this will be done by a separate function           & MANAGEMENT
      if false
         figure, hold on, box on, grid on
         plot(X(:),[Y1(:),Y2(:)])
         figure, hold on, box on, grid on
         plot([Y1(:),Y2(:)])
      else
         l_proc_trig(XRUN,YRUN,'run')
      end
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      % Trigger post-processing
      % do something with the data
      assignin('base','nitrgdata',bdatbuffer)
      disp('full trigger read and data sent to *nitrgdata* in base')
      % reset trigger text
      set(ud.handles{5},'string','')
      % reset trigger IDs
      ud.trig_cnt=0;  % reset trig counter
      ud.trig_idx=0;  % reset trig event index
   end
end


% GRAPHICS UPDATE
% ---------------
for iline=1:length(ud.handles{4})
   set(ud.handles{4}(iline), ...
      'xdata',timebuffer,'ydata',databuffer(:,iline))
end

% application update
ud.time=timebuffer;
ud.data=databuffer;
ud.bdat=bdatbuffer;
set(ud.handles{1},'userdata',ud)
return
% --------------------------------------------

function [udlis,hobj]=l_listbox(listname,hfig)
% listbox object to store individal test runs

% listbox data template
udlis=struct('counter',0, ...
   'dbid',[],'timestamp',{{}},'textnote',{{}},'dblock',{{}});
% listbox object
hobj=uicontrol('tag',listname,'parent',hfig,'style','listbox',...
   'unit','normalized','position',[0.92 0.11 0.075 0.815], ...
   'max',2,'userdata',udlis,'callback','');
% contextmenu control
hcm=uicontextmenu('parent',hfig,'userdata',hobj);
set(hobj,'uicontextmenu',hcm)
uimenu(hcm,'Label','Show sel','callback',{@l_proc_lst,'showsel'})
uimenu(hcm,'Label','Show sel avg','callback',{@l_proc_lst,'showselavg'})
uimenu(hcm,'Label','Export sel avg','callback',{@l_proc_lst,'expselavg'})
if exist('plot_cfitone','file')==2
   uimenu(hcm,'Label','Cfit one','callback',{@l_proc_lst,'cfitone'})
end
uimenu(hcm,'Label','Add note sel','callback',{@l_proc_lst,'addnote'})
uimenu(hcm,'Label','Delete sel','callback',{@l_proc_lst,'delsel'},'separator','on')
uimenu(hcm,'Label','Save UD ...','callback',{@l_proc_lst,'saveud'},'separator','on')
uimenu(hcm,'Label','Load UD ...','callback',{@l_proc_lst,'loadud'})
return
% --------------------------------------------

function l_window
% adjust window to current parameters

happ=findobj('tag','batch_ni_read_figure');
ud=get(happ,'userdata');

% acquisition and refresh parameters
Tfull=ud.dt_read;
Fevnt=ud.sa_evnt;
Fsamp=ud.sa_rate;
chans=ud.chan_id;
nblok=ud.nblocks;

% callback parameters (no safety checks)
lclickid=ud.fcn_lineclickfcnid;
fcnclick=ud.fcn_lineclickhans{lclickid};

% delete old lines
if ~isempty(ud.handles{4})
   delete(ud.handles{4})
end

% derived pars and preallocation
nsamp=Tfull*Fsamp;
if true % check maximum allowed size
   timebuffer=ones(nsamp,1)*NaN;
   databuffer=ones(nsamp,length(chans))*NaN;
else
   % do not allow large data objects
end

% NEW OBJECTS
set(ud.handles{2},'xlim',[0,Tfull]), zoom reset
hlin=[];
for ilin=1:length(chans)
   hlin(ilin)=line(ud.lspec{ilin}{:});
   set(hlin(ilin), ...
      'xdata',timebuffer,'ydata',databuffer(:,ilin), ...
      'buttondownfcn',{fcnclick});
end

% update application data
ud.time=timebuffer;
ud.data=databuffer;
ud.handles{4}=hlin;
set(ud.handles{1},'userdata',ud)
return
% --------------------------------------------

function l_init
% Initialize application

appname='NI signal reader';
apptag='batch_ni_read_figure';

hfig=figure('Name',appname,'NumberTitle','off','tag',apptag);
haxe=axes; set(haxe,'box','on','xgrid','on','ygrid','on'), hold on

% start/stop reading session
ht(1)=text('string','Start','units','normalized','position',[0.12,1.04], ...
   'buttondownfcn',{@batch_ni_read,'start'});
ht(2)=text('string','Stop','units','normalized','position',[0.19,1.04], ...
   'buttondownfcn',{@batch_ni_read,'stop'});
% init/prop/export/close session
ht(3)=text('string','InitNI','units','normalized','position',[0.52,1.04], ...
   'buttondownfcn',{@batch_ni_read,'initni'});
ht(4)=text('string','Poperties','units','normalized','position',[0.6,1.04], ...
   'buttondownfcn',{@batch_ni_read,'props'});
ht(5)=text('string','Export','units','normalized','position',[0.75,1.04], ...
   'buttondownfcn',{@batch_ni_read,'export'});
ht(6)=text('string','CloseNI','units','normalized','position',[0.87,1.04], ...
   'buttondownfcn',{@batch_ni_read,'closeni'});
% ... global UI props
set(ht,'fontsize',9,'backgroundcolor','w','edgecolor',[1,1,1]*0.35)
% ... trigger info
htt=text('string','abc','units','normalized','position',[0.99,0.99], ...
   'horizontalalignment','right','verticalalignment','top');

% application data
ud=l_init_ud(hfig,haxe,ht,htt);
set(hfig,'userdata',ud)
% --------------------------------------------

function ud=l_init_ud(hfig,haxe,ht,htt)
% Initialize data structure

Tspan=10;
Sa_rate=2048;     % other: 2048,3200,5120,6400
Nsa_event=256;
Nblocks=Sa_rate/Nsa_event;

ud=struct( ...
   ...    % NI objects
   'ni_device_obj',[], ...
   'ni_session_obj',[], ...
   'ni_listener_obj',[], ...
   ...    % basic reading parameters
   'sa_rate',Sa_rate, ...
   'dt_read',Tspan, ...
   'sa_evnt',Nsa_event, ...
   'axb_cal',{{[1,0]}}, ...
   ...    % state parameters
   'is_strt',false, ...
   'chan_id',[], ...
   'handles',{{hfig,haxe,ht,[],htt}}, ... % fig/axe/txt/lines/trig
   'lspec',{{{'color','b'},{'color','r'},{'color','m'},{'color','g'}}}, ...
   ...    % data read parameters
   't_end',0, ...
   't_abs',0, ...
   'sa_end',0, ...
   'bcount',0, ...
   'nblocks',Nblocks, ...
   'is_run',false, ...
   ...
   ...    % data containers
   'time',[], ...
   'data',[], ...
   ...
   ...    % triggers and callbacks
   'trig_flag',false, ...
   'trig_chan',0, ...
   'trig_lev',0.02, ...
   'trig_idx',0, ...
   'trig_cnt',0, ...
   'bdat',{cell(1,Nblocks)}, ...
   ...
   ...    % callback functions on events
   'fcn_lineclickfcnid',1, ...
   'fcn_lineclickhans',{{@l_proc_fft,@l_proc_frf}}, ...
   ...
   'fcn_evntcallfcnid',0, ...
   'fcn_evntcallhans',{{[],@l_proc_eswft}});
return
