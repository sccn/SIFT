function [fv2,loc] = downsamplesurf(fv,d,g)
if nargin<3, g = 0;end
N = size(fv.vertices,1);
d = round(d);
loc = 1:d:N;
fv2 = fv;
fv2.vertices = fv.vertices(loc,:);
dt = DelaunayTri(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3));
fv2.faces = freeBoundary(dt);
if g
    figure;patch('vertices',fv.vertices,'faces',fv.faces,'facecolor','g');title('Original Surface');
    figure;patch('vertices',fv2.vertices,'faces',fv2.faces,'facecolor','g');title('Reduced Surface');
end
