function guigr(varargin)
% GUIGR User interface for graph editing
%
% Usage: GUIGR to start a new graph
% GUIGR(NODES,EDGES) to continue editing a previous graph
% [No input checking implemented!]
%
% There are two editing modes: ADD and CHANGE
% to change between these modes, use the click on the first menu
% (either "add" or "change") or press Alt-A.
%
% ADD MODE: left click and drag to create new nodes and edges
% [snap to existing nodes]
%
% CHANGE MODE: allows deleting and shifting of nodes
% left click on a node to select (becomes marked with red dot)
% right click on selected node to delete or move
%
% To export Nodes and Edges to the base workspace click once onthe
% "export data" menu or press Alt-X.

% Created by H.R.A. Jagers, University of Twente, WL | Delft Hydraulics, The Netherlands
% on November 3rd, 1999

switch nargin,
case {0,2},
  TmpU=get(0,'units');
  set(0,'units','centimeter');
  Pos=get(0,'screensize');
  set(0,'units',TmpU);
  F=figure('doublebuffer','on','units','centimeter','position',[Pos(3)/2-5 Pos(4)/2-5 10 10],'units','pixels','resize','off');
  A=axes('parent',F,'units','normalized','position',[0 0 1 1],'xlim',[0 1],'ylim',[0 1],'dataaspectratio',[1 1 1],'box','on','xtick',[],'ytick',[],'visible','off');
  UD.Axes=A;
  UD.NodeH=line('xdata',[],'ydata',[],'linestyle','none','marker','.','parent',A,'clipping','off');
  UD.EdgeH=line('xdata',[],'ydata',[],'linestyle','-','parent',A,'clipping','off');
  UD.Cord=line('xdata',[],'ydata',[],'linestyle','-','parent',A,'erasemode','xor','color',[.5 .5 .5],'clipping','off');
  UD.SelH=line('xdata',[],'ydata',[],'marker','.','markersize',18,'color','r','parent',A,'visible','off','clipping','off');
  UD.Selected=0;
  UD.StartCord=0;
  UD.Nodes=[];
  UD.Edges=[];
  UD.Mode='add';
  uimenu('parent',F,'label','ch&ange','callback','guigr change');
  uimenu('parent',F,'label','e&xport data','callback','guigr export');
  UD.EditMenu=uicontextmenu;
  set(UD.SelH,'uicontextmenu',UD.EditMenu);
  uimenu(UD.EditMenu,'label','&move','callback','guigr move');
  uimenu(UD.EditMenu,'separator','on','label','&delete','callback','guigrdelete');
  uimenu(UD.EditMenu,'label','delete &edge','callback','guigrdeledge','enable','off');
  if nargin==2, % draw graph
    UD.Nodes=varargin{1};
    set(UD.NodeH,'xdata',UD.Nodes(:,1),'ydata',UD.Nodes(:,2));
    UD.Edges=varargin{2};
    I=transpose([ones(size(UD.Edges,1),1) UD.Edges]);
    I=UD.Nodes(I(:),:);
    I(1:3:end,:)=NaN;
    set(UD.EdgeH,'xdata',I(:,1),'ydata',I(:,2));
  end;
  set(F,'windowbuttondownfcn','guigrdown','userdata',UD,'menubar','none','name','[add mode]','numbertitle','off');
case 1,
  Pnt=get(getudf(gcbf,'Axes'),'currentpoint');
  Pnt=Pnt(1,1:2);
  UD=get(gcbf,'userdata');
  switch varargin{1},
  case 'motion',
    Pos=get(gcbf,'position');
    T=get(0,'pointerlocation');
    T=min(Pos(1:2)+Pos(3:4),max(Pos(1:2),T));
    set(0,'pointerlocation',T);
    Pnt=get(getudf(gcbf,'Axes'),'currentpoint');
    Pnt=Pnt(1,1:2);
    DistSq=(UD.Nodes(:,1)-Pnt(1)).^2+(UD.Nodes(:,2)-Pnt(2)).^2;
    [minDistSq,i]=min(DistSq);
    if minDistSq<0.001,
      Pnt=UD.Nodes(i,:);
    end;
    switch UD.Mode,
    case 'add',
      X=get(UD.Cord,'xdata'); X(2)=Pnt(1);
      Y=get(UD.Cord,'ydata'); Y(2)=Pnt(2);
      set(UD.Cord,'xdata',X,'ydata',Y);
    case 'change',
      X=get(UD.Cord,'xdata'); X(2:3:end)=Pnt(1);
      Y=get(UD.Cord,'ydata'); Y(2:3:end)=Pnt(2);
      set(UD.Cord,'xdata',X,'ydata',Y);
      set(UD.SelH,'xdata',Pnt(1),'ydata',Pnt(2));
    end;
  case 'change',
    if strcmp(get(gcbf,'windowbuttonupfcn'),''),
      set(gcbo,'label','&add','callback','guigr add');
      set(gcbf,'windowbuttondownfcn','guigr select','name','[change mode]');
      UD.Selected=0;
      UD.Mode='change';
    end;
  case 'add',
    set(gcbo,'label','ch&ange','callback','guigr change');
    set(gcbf,'windowbuttondownfcn','guigr down','name','[add mode]');
    set(UD.SelH,'visible','off');
    UD.Mode='add';
  case 'export',
    assignin('base','Nodes',UD.Nodes);
    assignin('base','Edges',sort(UD.Edges,2));
  case 'select',
    switch get(gcbf,'selectiontype'),
    case 'normal',
      if isempty(UD.Nodes),
        minDistSq=100;
      else,
        DistSq=(UD.Nodes(:,1)-Pnt(1)).^2+(UD.Nodes(:,2)-Pnt(2)).^2;
        [minDistSq,i]=min(DistSq);
      end;
      if minDistSq<0.001,
        Pnt=UD.Nodes(i,:);
        set(UD.SelH,'visible','on','xdata',Pnt(1),'ydata',Pnt(2));
        UD.Selected=i;
      end;
    end;
  case 'deledge',
    UD.StartCord=UD.Selected;
    set(UD.Cord,'xdata',[UD.Nodes(UD.Selected,1) Pnt (1)],'ydata',[UD.Nodes(UD.Selected,2) Pnt(2)],'visible','on');
  case 'delete',
    UD.Nodes(UD.Selected,:)=[];
    set(UD.NodeH,'xdata',UD.Nodes(:,1),'ydata',UD.Nodes(:,2));
    set(UD.SelH,'visible','off');
    UD.Edges(any(UD.Edges==UD.Selected,2),:)=[];
    if isempty(UD.Edges),
      set(UD.EdgeH,'xdata',[],'ydata',[]);
    else,
      I=UD.Edges>UD.Selected;
      UD.Edges(I)=UD.Edges(I)-1;
      I=transpose([ones(size(UD.Edges,1),1) UD.Edges]);
      I=UD.Nodes(I(:),:);
      I(1:3:end,:)=NaN;
      set(UD.EdgeH,'xdata',I(:,1),'ydata',I(:,2));
    end;
    UD.Selected==0;
  case 'move',
    if isempty(UD.Edges),
      I=[];
    else,
      I=any(UD.Edges==UD.Selected,2);
    end;
    UD.StartCord=UD.Selected;
    if any(I),
      I=UD.Edges(I,:);
      I=I(I(:)~=UD.Selected);
      I=transpose([ones(size(I,1),2) I]);
      I(2,:)=UD.Selected;
      I=UD.Nodes(I(:),:);
      I(1:3:end,:)=NaN;
      set(UD.Cord,'xdata',I(:,1),'ydata',I(:,2),'visible','on');
    else,
      set(UD.Cord,'xdata',[],'ydata',[],'visible','off');
    end;
    set(UD.SelH,'erasemode','xor');
    set(gcbf,'windowbuttonmotionfcn','guigr motion', ...
             'windowbuttonupfcn','guigr up', ...
             'windowbuttondownfcn','');
  case 'down',
    switch get(gcbf,'selectiontype'),
    case 'normal',
      if isempty(UD.Nodes),
        minDistSq=100;
      else,
        DistSq=(UD.Nodes(:,1)-Pnt(1)).^2+(UD.Nodes(:,2)-Pnt(2)).^2;
        [minDistSq,i]=min(DistSq);
      end;
      if minDistSq<0.001,
        Pnt=UD.Nodes(i,:);
      else,
        i=size(UD.Nodes,1)+1;
        UD.Nodes(i,1:2)=Pnt;
        set(UD.NodeH,'xdata',UD.Nodes(:,1),'ydata',UD.Nodes(:,2));
      end;
      UD.StartCord=i;
      set(gcbf,'windowbuttonmotionfcn','guigr motion', ...
               'windowbuttonupfcn','guigr up', ...
               'windowbuttondownfcn','');
      set(UD.Cord,'xdata',[Pnt(1) Pnt(1)],'ydata',[Pnt(2) Pnt(2)],'visible','on');
    case 'alt',
    end;
  case 'up',
    Pos=get(gcbf,'position');
    T=get(0,'pointerlocation');
    T=min(Pos(1:2)+Pos(3:4),max(Pos(1:2),T));
    set(0,'pointerlocation',T);
    Pnt=get(getudf(gcbf,'Axes'),'currentpoint');
    Pnt=Pnt(1,1:2);
    switch UD.Mode,
    case 'add',
      set(gcbf,'windowbuttonmotionfcn','', ...
               'windowbuttonupfcn','', ...
               'windowbuttondownfcn','guigr down');
      set(UD.Cord,'visible','off');
      DistSq=(UD.Nodes(:,1)-Pnt(1)).^2+(UD.Nodes(:,2)-Pnt(2)).^2;
      [minDistSq,i]=min(DistSq);
      if minDistSq<0.001,
        Pnt=UD.Nodes(i,:);
      else,
        i=size(UD.Nodes,1)+1;
        UD.Nodes(i,1:2)=Pnt;
        set(UD.NodeH,'xdata',UD.Nodes(:,1),'ydata',UD.Nodes(:,2));
      end;
      if UD.StartCord~=i,
        Drawn=ismember_bc([UD.StartCord i; i UD.StartCord],UD.Edges,'rows');
        if ~any(Drawn),
          Ex=get(UD.EdgeH,'xdata'); Ex(size(Ex,2)+(1:3))=[NaN UD.Nodes(UD.StartCord,1) Pnt(1)];
          Ey=get(UD.EdgeH,'ydata'); Ey(size(Ey,2)+(1:3))=[NaN UD.Nodes(UD.StartCord,2) Pnt(2)];
          UD.Edges(size(UD.Edges,1)+1,1:2)=[i UD.StartCord];
          UD.StartCord=0;
          set(UD.EdgeH,'xdata',Ex,'ydata',Ey);
        end;
      end;
    case 'change',
      set(gcbf,'windowbuttonmotionfcn','', ...
               'windowbuttonupfcn','', ...
               'windowbuttondownfcn','guigr select');
      set(UD.Cord,'visible','off');
      DistSq=(UD.Nodes(:,1)-Pnt(1)).^2+(UD.Nodes(:,2)-Pnt(2)).^2;
      [minDistSq,i]=min(DistSq);
      if minDistSq<0.001,
        UD.Nodes(UD.Selected,:)=[];
        set(UD.NodeH,'xdata',UD.Nodes(:,1),'ydata',UD.Nodes(:,2));
        UD.Edges(UD.Edges==UD.Selected)=i;
        UD.Edges(UD.Edges(:,1)==UD.Edges(:,2),:)=[];
        if isempty(UD.Edges),
          set(UD.EdgeH,'xdata',[],'ydata',[]);
        else,
          I=UD.Edges>UD.Selected;
          UD.Edges(I)=UD.Edges(I)-1;
          I=transpose([ones(size(UD.Edges,1),1) UD.Edges]);
          I=UD.Nodes(I(:),:);
          I(1:3:end,:)=NaN;
          set(UD.EdgeH,'xdata',I(:,1),'ydata',I(:,2));
        end;
        if i>UD.Selected,
          UD.Selected=i-1;
        else,
          UD.Selected=i;
        end;
      else,
        i=UD.StartCord;
        UD.Nodes(i,:)=Pnt;
        set(UD.NodeH,'xdata',UD.Nodes(:,1),'ydata',UD.Nodes(:,2));
        set(UD.SelH,'xdata',Pnt(1),'ydata',Pnt(2));
        I=transpose([ones(size(UD.Edges,1),1) UD.Edges]);
        I=UD.Nodes(I(:),:);
        I(1:3:end,:)=NaN;
        set(UD.EdgeH,'xdata',I(:,1),'ydata',I(:,2));
      end;
      set(UD.SelH,'erasemode','normal');
    end;
  end;
  set(gcbf,'userdata',UD);
end;