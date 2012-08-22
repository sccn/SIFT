function figure_wire(surf,Fcolor);
% it colors mesh faces and eges using Fcolor and Ecolor

figure
Ecolor=[0.5 0.5 0.5]; %color of wire frame

p=patch(surf);
set(p,'FaceColor',Fcolor,'EdgeColor',Ecolor);
daspect([1 1 1]);
view(3); axis tight
camlight; 
lighting gouraud
alpha(0.9);
