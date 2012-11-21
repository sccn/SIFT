function ICAshow(sigMatrix, varargin);
% ICASHOW ---- plot the input signals with optional title and comments
%
% ���ø�ʽ��
%    ICAshow(sigMatrix,'title','title of picture','comment','comment of the picture');
%            �����źŵĲ��Σ������Ϊ'title of picture'��ͼ���µ�˵��Ϊ'comment of the picture'
%
% ����˵�� ......
%      'title'  ---- ͼƬ�ı���
%    'comment'  ---- ��ͼƬ��˵��
%
% ���ߣ������֣�Zhang Zhi-Lin��
%       zlzhang@uestc.edu.cn
% �汾��1.0
% ���ڣ�2003��11��1��


flagcomment=0;   % ��־�ſ�ѡ����'comment'�Ƿ�ѡ�����Ϊ1��ʾ�û�ѡ��ò���������ò���ֵ 
flagtitle=0;     % ��־�ſ�ѡ����'title'�Ƿ�ѡ�����Ϊ1��ʾ�û�ѡ��ò���������ò���ֵ 
% ��ȡ��ѡ����
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs.\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'title'
                titlechar=varargin{i+1};
                flagtitle=1;
            case 'comment'
                commentchar=varargin{i+1};
                flagcomment=1;
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

rows=size(sigMatrix,1);
cols=size(sigMatrix,2);
%sigMatrix = remstd(sigMatrix);

figure;
for i = rows : -1 : 1
    subplot(rows,1,i);
    plot(sigMatrix(i,:)); axis tight;
    
    if  flagcomment==1 & i == rows
        xlabeltext=['\bf',commentchar];
        xlabel(xlabeltext);
    end
    
    if  flagtitle==1 & i ==1
        titletext=['\bf',titlechar];
        title(titletext);
    end
     
     %axis([-inf,inf,-5,5]);
end


function sig = remstd(OldSig)
% remstd --- make the input signals' variance unit
% The input OldSig is in the form N*P, where N is the number of signals, and P is
% the length of each signal, then the returned sig is also in the form N*P,
% but each row (each signal) has unit variance.
%
% Command:
%    sit = remstd( OldSig );
%
% See also:
%     centering    whiten    standarize
%
% Author��Zhi-Lin Zhang
%         zlzhang@uestc.edu.cn
%         http://teacher.uestc.edu.cn/teacher/teacher.jsp?TID=zzl80320
%
% Last version: 1.5     Date: Apr.17,2005
%      Version: 1.0     Date: Oct.31,2003


n=size(OldSig,1);
for t=1:n
    sig(t,:)=OldSig(t,:)/std(OldSig(t,:));
end



