%%
run_readlsl('mystream')

%% BCILAB streamviewer

f=figure; %('KeyPressFcn',@(varargin) evalin('base',strcmp(varargin{2}.Key,'uparrow'));

spacing = 30;
timerange = 10;
while ~isempty(get(f)) 
    tmp=onl_peek('mystream',timerange);
    plot(bsxfun(@plus,tmp.data',(0:tmp.nbchan-1)*spacing-mean(tmp.data'))); 
    axis([0 tmp.pnts -300 tmp.nbchan*spacing+300]);
    set(gca,'YTick',(0:tmp.nbchan-1)*spacing,'YTickLabel',{tmp.chanlocs.labels});
    drawnow;
end

%% Barplot
display_range = 262000;
tmp=onl_peek('mystream',1); 
f=figure;
while ~isempty(get(f)) 
    tmp=onl_peek('mystream',1); 
    bar(1:tmp.nbchan, mean(tmp.data')/256,'Horizontal','on'); 
    axis([-display_range display_range 1 tmp.nbchan]);
    %bar(1:tmp.nbchan, log(std(tmp.data')),'Horizontal','on'); 
    set(gca,'YTick',1:tmp.nbchan,'YTickLabel',{tmp.chanlocs.labels});
    drawnow; 
end

%%
figure;
set(gcf,'Colormap',rand(3,3));
bar(1:10,rand(10,3));