function elem=removedupelem(elem)
%
% elem=removedupelem(elem)
%
% remove doubly duplicated elements
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
%
% input:
%    elem: list of elements (node indices)
%
% output:
%    elem: element list after removing the duplicated elements
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

[el,count1,count2]=unique_bc(sort(elem')','rows');
bins=hist(count2,1:size(elem,1));
cc=bins(count2);
elem(find(cc>0&mod(cc,2)==0),:)=[];
