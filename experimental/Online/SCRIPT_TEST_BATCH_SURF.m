
%% render each cylinder through a separate call to surf in a loop

figloop=figure;
hax=axes('units','normalized','drawmode','fast');
axis([0 1 0 1]);

K = 100;

x = rand(2,K);
y = rand(2,K);
z = rand(2,K);
colors = rand(3,K);

tic 
hold on
for i=1:K
    [xc yc zc] = cylinder2P(0.01, 11, 2, [x(1,i) y(1,i) z(1,i)], [x(2,i) y(2,i) z(2,i)]);
    colorarray = repmat(reshape(colors(:,i), 1,1,3), [size(zc,1) size(zc,2) 1]);
    handles = surf(xc, yc, zc, colorarray,'edgecolor', 'none', ...
        'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', 'phong');
end
hold off
camlight headlight
drawnow;

t=toc;

axis([0 1 0 1]);
fprintf('loop time = %0.5g\n',t);

%% render all cylinders at once

clear xc zc yc colorarray;

figbatch=figure;
hax=axes('units','normalized','drawmode','fast');


[xc yc zc] = deal(zeros(3*K,11));
colorarray = zeros(size(xc,1), 11, 3);

tic 
q = 1;
for i=1:K
    % make the cylinder
    [xc(q:q+1,:) yc(q:q+1,:) zc(q:q+1,:)] = cylinder2P(0.01, 11, 2, [x(1,i) y(1,i) z(1,i)], [x(2,i) y(2,i) z(2,i)]);
    colorarray(q:q+2,:,:) = repmat(reshape(colors(:,i), 1,1,3), [3 11 1]);
    q = q+2;
    % insert a row of NaNs
    [xc(q,:) yc(q,:) zc(q,:)] = deal(nan(1,size(xc,2)));
    q = q+1;
end

% render all cylinders at once
handles = surf(xc, yc, zc,colorarray,'edgecolor', 'none', ...
    'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', 'phong');
camlight headlight
drawnow;

t=toc;

axis([0 1 0 1]);
fprintf('batch time = %0.5g\n',t);


% [xc yc zc] = cylinder2P(3, 11, 2, [x(1,i) y(1,i) z(1,i)], [x(2,i) y(2,i) z(2,i)]);
% colorarray = repmat(reshape([155 10 40], 1,1,3), [size(zc,1) size(zc,2) 1]);
% handles = surf(hax,xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
%     'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', 'phong');