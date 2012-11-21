function [S,W,E] = FastICA(Z,varargin)
% FastICA ---- Fast Fixed-Point ICA algorithm
%     ������İ׻����źţ����ֵ��λ����ͣ��Լ���ز���������symmetric����������
%     deflation������������ָ���ԭʼ�ź�S�Լ���Ӧ�ļ�Ȩ����W��������symmetric����ʱ��
%     ������ٸ��۲��ź�����Ƴ����ٸ�����������������deflation����ʱ�����Թ���������Ŀ�Ķ�������
%     ������Ķ���������Ŀ��С�ڹ۲��ź���Ŀ��
%  
% ����:
%     [S,W] = FastICA(Z,'method','defl',  'origW',W0,  'NonlinFun','tanh',  'a1',1,  'a2',1, 
%             'maxNumIteration',300,  'OverValue',0.001,  'report','off',  'mixMatrix',VA,  'ICNum',3)
%
% ����˵��......
%       ��������    ||   ����ȡֵ       || ����˵��                   
%               Z                       : �׻��˵��źţ�ÿ�д���һ���ź�
%         'method'    ['defl'|'symm']   : ���õ�����������
%                                         'defl'Ϊdeflationary����
%                                         'symm'Ϊsymmetric������Ĭ��ֵΪ'symm'                                                mm'
%          'origW'                      : �����ʼ�ļ�Ȩ����
%                                         ��ʡ�ԣ�����ϵͳ������ʼ����������
%      'NonlinFunc'    ['pow3'|         : ���õķ����Ժ�����'pow3'Ϊg(y)=y^3
%                       'tanh'|           'tanh'Ϊg(y)=tanh(a1*y)  
%                       'gaus']           'gaus'Ϊg(y)=y*exp(-y^2/2)                                                            
%              'a1'                     : �����Ժ��� 'tanh'�е��õ�ֵ,1<=a1<=2��Ĭ��ֵ1
%              'a2'                     : �����Ժ��� 'gaus'�е��õ�ֵ��Ĭ��ֵ1
% 'maxNumIteration'                     : ���ĵ���������Ĭ��ֵΪ500
%       'OverValue'                     : ����������������ֵ����ԽС������W*W'Խ���ڵ�λ��
%                                         Ĭ��ֵΪ0.0001
%          'report'    ['on'|'off']     : �Ƿ���Ҫ��ӡ������������Ϣ
%                                         'on'Ϊ��ӡ��'off'Ϊ����ӡ��Ĭ��Ϊ'on'
%       'mixMatrix'                     : �� V*A,VΪ�׻�����AΪ��Ͼ���
%           'ICNum'                     : ����Ĺ����źŸ�����ע��˲�������'method'����Ϊ'defl'ʱʹ��
%                E                      : ��¼ÿ�ε����������ָ��������ָ���ɺ���CTPI(W'*V*A,'norm')����
%                S                      : �ָ����ź�: S=W'* Z
%                W                      : �ָ��źŵļ�Ȩ����,W=[b1,b2,...bm]
%
% ���ø�ʽ��
%     �� [S,W]=FastICA(WhitenedSig,'nonlinfunc','POW3','maxNumIteration',100) 
%     ����Ŀ�ѡ��������ɶԡ����ǲ������������Ƕ�Ӧ�Ĳ���ֵ��ע����������ֵ�����ַ�����
%
% ���ߣ������֣�Zhang Zhi-Lin��
%       zlzhang@uestc.edu.cn          zzl.private@eyou.com
%       http://teacher.uestc.edu.cn/teacher/teacher.jsp?TID=zzl80320
% 
% �汾��1.0    ���ڣ�2003��10��31��
% �汾��1.1    ���£�2004��12��20��
% �汾��1.2    ���£�2005�� 7�� 1��
% �汾��1.3    ���£�2006�� 4��11��   ��������deflation����ʱ���ɹ��Ƶ��ź���Ŀ���ص���Դ�ź���Ŀ��
%
% �ο����ף�
% [1] A.Hyvarinen. Fast and Robust Fixed-Point Algorithms for Independent
%     Component Analysis. IEEE Trans. on Neural Networks,10(3):626-634,1999
%




% ===================================================================
% ȡĬ��ֵ
ICNum=size(Z,1);    % ����Ϊ���ٸ��źţ���ָ������ٸ��ź�    
method='symm';
NonlinFunc='tanh';
a1=1;
a2=1;
maxNumIteration=500;
mixMatrix = eye(ICNum);
OverValue=0.0001;
report='on';
numSamples=size(Z,2);    % ÿ���źŵĲ������������źŵĳ���

UserInputW=0;   % �ж��Ƿ����û������ʼ�ļ�Ȩ����Ϊ1�������û����롣

% ===================================================================
% ��ȡ��ѡ����
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'method'
                method=varargin{i+1}; method=lower(method);
            case 'origw'  % origW
                W=varargin{i+1}; 
                UserInputW=1;   % ���û������ʼ�ļ�Ȩ����
            case 'nonlinfunc'   % NonlinFunc
                NonlinFunc=varargin{i+1}; NonlinFunc=lower(NonlinFunc);
            case 'a1'
                a1=varargin{i+1};  
            case 'a2'
                a2=varargin{i+1};  
            case 'maxnumiteration'   % maxNumIteration
                maxNumIteration=varargin{i+1}; 
            case 'mixmatrix'
                mixMatrix = varargin{i+1};
            case 'overvalue'         % OverValue
                OverValue=varargin{i+1};  
            case 'icnum'             % ICNum
                ICNum = varargin{i+1};
            case 'report'
                report=varargin{i+1};  report=lower(report);
            otherwise
                % ������������д���
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

if UserInputW == 0
    % ��ϵͳ������ʼ��������Ȩ����W
    %W=orth(rand(ICNum)-0.5);
    W=orth(eye(ICNum));
end

% ===================================================================
% ��������������ȷ��
if ~( strcmp(report,'on') | strcmp(report,'off') )
    error(sprintf('Illegal value [ %s ] for parameter: ''report''\n', report));
end

if ~isreal(Z)
    error('Input matrix has an imaginary part.');
end

if ~( strcmp(method,'symm') | strcmp(method,'defl'))
    error(sprintf('Illegal value [ %s ] for parameter: ''method''\n', method));
end

if ~( strcmp(NonlinFunc,'pow3') | strcmp(NonlinFunc,'tanh')| ...
        strcmp(NonlinFunc,'gaus') | strcmp(NonlinFunc,'skew') )
    error(sprintf('Illegal value [ %s ] for parameter: ''NonlinFunc''\n', NonlinFunc));
end

if ~isnumeric(a1)
    error(sprintf('Illegal value for parameter: ''a1'',it must be number\n', a1));
end

if ~isnumeric(a2) 
    error(sprintf('Illegal value for parameter: ''a2'',it must be number\n', a2));
end

if (~isnumeric(maxNumIteration) ) |  maxNumIteration < 1 
    error(sprintf('Illegal value for parameter: ''maxNumIteration''\n', maxNumIteration));
end

if ~isnumeric(OverValue) 
    error(sprintf('Illegal value for parameter: ''OverValue'', it must be number.\n', OverValue));
end

if ( size(W,1)~=size(W,2) ) | ( ~isreal(W) ) | size(W,1)~=ICNum 
    error(sprintf('Illegal value for parameter: ''origW''\n'));
end

% ===================================================================
% ��ӡ������Ϣ
if strcmp(report,'on') 
    fprintf('\n=======================================\n');
    fprintf('   Information about parameters...\n');
    fprintf('=======================================\n');
    fprintf('Number of input signals: %d\n',size(Z,1));
    fprintf('Number of estimated signals: %d\n',ICNum);
    fprintf('Length of signals:%d\n',numSamples);    
    fprintf('method = %s\n',method);
%     fprintf('initial W = \n');disp(W);
    fprintf('NonlinFunc = %s\n',NonlinFunc);
    fprintf('a1 = %d\n',a1);
    fprintf('a2 = %d\n',a2);
    fprintf('maxNumIteration = %d\n',maxNumIteration);
    fprintf('OverValue = %d\n',OverValue);
    fprintf('report = %s\n',report);
%     fprintf('mixMatrix = \n'); disp(mixMatrix);
end

if ICNum<size(Z,1) & strcmp(method,'symm')
    fprintf('$$ Warning: You use the ''symm'' method, but the number of estimated signals is less than that of the observed signals ! $$\n');
end
 
if strcmp(report,'on') 
    fprintf('\n=======================================\n');
    fprintf('          Starting FastICA  ...\n'); 
    fprintf('=======================================\n');
end

% ===================================================================
% symmetric ��������
if strcmp(method,'symm')
    if( strcmp(report,'on') )   fprintf('Using symmetric method \n'); end
    
    WOld=zeros(size(W));  % WOld����������һ�ε�W��ֵ���Ա�����ڵ�Wֵ��Ƚ�
    
    %--------------------------------------------
    % fast fixed-point ICA �㷨�ĵ���
    for loop=1:maxNumIteration + 1
        
        % ��������������maxNumIteration����û�����������ӡ��Ϣ�����ص�ǰ��W��Sֵ��
        if loop == maxNumIteration + 1  
            fprintf('Cannot convergence even after %d steps \n',maxNumIteration);
            S=W'*Z;
            return;
        end
        
        % Symmetric ������
        W=W*real(inv(W'*W)^(1/2)); 
        
        % ��¼���ε���������ָ��
        E(loop) = CTPI(W'* mixMatrix,'norm');
        
        % �ж��Ƿ���ϵݹ�����������������Ҳ���ǵ��˷��������෴������
        minAbsCos = min(abs(diag(W'*WOld)));  
        
        if(strcmp(report,'on'))
            % ��ӡ�������̵���Ϣ
            fprintf('-------------------\n');
            fprintf('Step No.%d, change in value of estimate: %.4f \n',loop,1-minAbsCos);           
        end
        
        if( 1-minAbsCos < OverValue)
            if( strcmp(report,'on') )  
                fprintf('\n=======================================\n');
                fprintf('  Convergence after %d steps\n',loop);
                fprintf('=======================================\n');
            end
            S=W'*Z;
            return;                %��������ѭ������������������
        end
        
        WOld=W;
        
        % ���ò�ͬ�ķ����Ժ�������
        switch NonlinFunc
            % pow3
            case 'pow3'
                W=(Z*((Z'*W).^3))/numSamples-3*W;   % �������۲���������ƽ��
            case 'tanh'
                hypTan=tanh(a1*Z'*W);
                W=Z*hypTan/numSamples -ones(size(W,1),1)*sum(1-hypTan.^2).*W/numSamples*a1;       
            case 'gaus'
                Y=Z'*W;
                Y2=Y.^2;
                ex=exp(-a2*Y2/2);
                gauss=Y.*ex;
                degauss=(1-a2*Y2).*ex;
                W=Z*gauss/numSamples - ones(size(W,1),1)*sum(degauss).*W/numSamples;
        end
        
        
    end   % ��������
end   % symmetric method ����


% ===================================================================
% deflation ��������,��Ȩ����B��ÿ��������һ����һ���������
if strcmp(method,'defl')
    if( strcmp(report,'on') )   fprintf('Using deflation method \n'); end
    
    B = zeros(size(Z,1),ICNum);   % ���Ҫ����ļ�Ȩ����
    
    % �����ʼ�ļ�Ȩ����
    if UserInputW==0  % ��ϵͳ������ʼ�ļ�Ȩ����
        w=rand(size(Z,1),1)-0.5;
    else    % ���û������ʼ�ļ�Ȩ����
        w=W(:,1);
    end
    w=w-B*B'*w;     % ������
    w=w/norm(w);    % ��һ��
    wOld=zeros(size(w));   % wOld����������һ�ε�wֵ
    
    %--------------------------------------------
    for i=1:ICNum                % �������ICNum����������
        
        % ------------------------------
        % fast fixed-point ICA �㷨�ĵ���
        for loop=1:maxNumIteration+1
            % �Ѿ��ҵ��Ļ�����������һ�������ռ䡣���ڰѵ�ǰ������ͶӰ����������ռ��
            % �����ռ���Ӷ�ʵ����������ע�⵽֮�����ܹ����þ���B��ȷ����ͶӰ������
            % Ϊ����B�е���������ͶӰ�������κ�Ӱ�졣
            w=w-B*B'*w;     % ������
            w=w/norm(w);    % ��һ��
            
            if loop == maxNumIteration+1     % ����������maxNumIterationʱ�����˳�����
                fprintf('IC No.%d cannot converge in %d iteration.\n',i,maxNumIteration);
                S=B'*Z;
                return;
            end
            
            % �ж��Ƿ�������������ʱ��w��wOld������ͬ�����������෴��
            % ���ԣ��ж����������������������
            if norm(w-wOld) < OverValue | norm(w+wOld) < OverValue
                if( strcmp(report,'on') )  
                    fprintf('IC No.%d  ---------- Computed after %d iteration\n',i,loop-1);
                end
                B(:,i)=w;  % ����ǰ��wֵ������B��
                break;
            end
            
            wOld=w;
            
            % ���ò�ͬ�ķ����Ժ�������
            switch NonlinFunc
                
                case 'pow3'
                    w = ( Z*( (Z'*w).^3 ) )/numSamples - 3*w;
                    
                case 'tanh'
                    hypTan = tanh(a1*Z'*w);
                    w = (Z*hypTan - a1*sum(1-hypTan.^2)'*w)/numSamples;
                    
                case 'gaus'
                    u=Z'*w;
                    u2=u.^2;
                    ex=exp(-a2*u2/2);
                    gauss=u.*ex;
                    degauss=(1-a2*u2).*ex;
                    w=(Z*gauss - sum(degauss)'*w)/numSamples;
            end
         
        end  % ��B��ĳ��������wѭ��maxNumIteration �� 
       
    end  % ��ICNum��������������
    
    W=B;     % ��Ȩ��������������W
    S=W'*Z;  % �ָ����ź�
    
    if ICNum == size(Z,1)
        E = CTPI(W'* mixMatrix,'norm');   % ����ָ��
    end
    
end % deflation method ����



%------------------------------------------------------------------------
function E = CTPI(P,varargin)
% CTPI -- Returns the cross-talk error, which is used in blind source
%         separation to measure the accuracy of separation. The larger it
%         is, the worse performance of separation.
%
% Command:
%          E = CTPI(P,'sum');    
%          E = CTPI(P,'ave');              
%          E = CTPI(P,'norm'); 
%          E = CTPI(P,'single');
%
% Parameter:
% < Input >
%            P --- P=WA,where W is demixed matrix and A is mixing matrix
%                  and whose P(i,j) is the ijth element of the N*N matrix P
%        'sum' --- compute the performance index as follows:
%                          N      N     |P(i,j)|               N      N     |P(j,i)|         
%              CTPI(P) =  sum (  sum --------------- - 1  ) + sum (  sum --------------- - 1  )
%                         i=1    j=1  max(|P(i,j)|)           j=1    i=1  max(|P(j,i)|)
%
%        'ave' --- compute the performance index as follows:
%                         1   N      N     |P(i,j)|              1   N      N     |P(j,i)|         
%              CTPI(P) = --- sum (  sum --------------- - 1 ) + --- sum (  sum --------------- - 1 )
%                         N   i=1   j=1  max(|P(i,j)|)           N  j=1    i=1  max(|P(j,i)|)
%
%     'single' --- compute the performance index as follows:
%                           1   N      N     |P(i,j)|
%              CTPI(P) =   --- sum (  sum --------------- - 1  )
%                           N  i=1    j=1  max(|P(i,j)|)
%
%       'norm' --- compute the performance index as follows:
%                                               max(|P(i,j)|^2)    max(|P(j,i)|^2) 
%                            1        1   m    1<=j<=m            1<=j<=m
%              CTPI(P) =   ----- (m- --- sum ( ---------------- + -----------------) )
%                           m-1       2  i=1     m                 m                     
%                                               sum(|P(i,j)|^2)   sum(|P(j,i)|^2)
%                                               j=1               j=1
%                  the performance metric has the following features:
%                  (1) lies in [0,1] for all matrices P
%                  (2) E = 1 if and only if |Pij|^2 = |Ppq|^2 for all i,j,p,q in
%                      the range [1,m] (i.e.,maximally mixed sources in the system
%                      outputs)
%                  (3) E = 0 if and only if P = \Phi * D (i.e., separated sources
%                      in the system outputs)
%
% Author  : Zhi-Lin Zhang ( zlzhang@uestc.edu.cn )
% Version : 1.0
% Date    : June 13, 2005 



[N,J] = size(P);
if N ~= J
    fprintf('Number of columns are NOT equal to that of rows!\nExit...\n');
    return;
end

if strcmp(lower(varargin{1}),'norm')  
    % compute the performance:
    %                                     max(|P(i,j)|^2)    max(|P(j,i)|^2) 
    %                  1        1   N    1<=j<=N            1<=j<=N
    %    CTPI(P) =   ----- (N- --- sum ( ---------------- + -----------------) )
    %                 N-1       2  i=1     N                 N                     
    %                                     sum(|P(i,j)|^2)   sum(|P(j,i)|^2)
    %                                     j=1               j=1
    P = P.^2;
    E =  (N-sum( max(P)./sum(P) + max(P')./sum(P') )/2 ) / (N-1)  ;
    
elseif strcmp(lower(varargin{1}),'sum')  
    %    'sum' --- compute the performance index as follows:
    %                      N      N     |P(i,j)|               N      N     |P(j,i)|         
    %          CTPI(P) =  sum (  sum --------------- - 1  ) + sum (  sum --------------- - 1  )
    %                     i=1    j=1  max(|P(i,j)|)           j=1    i=1  max(|P(j,i)|)   
    P = abs(P);
    sum1 = 0;
    for i = 1:N
        maxEle = max(P(i,:));
        sum1 = sum1 + sum(P(i,:))/maxEle;
    end
    sum1 = sum1 - N;
    
    sum2 = 0;
    for j = 1:N
        maxEle = max(P(:,j));
        sum2 = sum2 + sum(P(:,j))/maxEle;
    end
    sum2 = sum2 - N;
    
    E = sum1 + sum2;
    
elseif strcmp(lower(varargin{1}),'ave')  
%    'ave' --- compute the performance index as follows:
%                 1   N      N     |P(i,j)|              1   N      N     |P(j,i)|         
%      CTPI(P) = --- sum (  sum --------------- - 1 ) + --- sum (  sum --------------- - 1 )
%                 N   i=1   j=1  max(|P(i,j)|)           N  j=1    i=1  max(|P(j,i)|)
%
    P = abs(P);
    sum1 = 0;
    for i = 1:N
        maxEle = max(P(i,:));
        sum1 = sum1 + sum(P(i,:))/maxEle;
    end
    sum1 = sum1 - N;
    
    sum2 = 0;
    for j = 1:N
        maxEle = max(P(:,j));
        sum2 = sum2 + sum(P(:,j))/maxEle;
    end
    sum2 = sum2 - N;
    
    E = (sum1 + sum2)/N;
    
elseif strcmp(lower(varargin{1}),'single')  
%     'single' --- compute the performance index as follows:
%                           1   N      N     |P(i,j)|
%              CTPI(P) =   --- sum (  sum --------------- - 1  )
%                           N  i=1    j=1  max(|P(i,j)|)
    P = abs(P);
    sum1 = 0;
    for i = 1:N
        maxEle = max(P(i,:));
        sum1 = sum1 + sum(P(i,:))/maxEle;
    end
    E = (sum1 - N)/N;
    
else
    fprintf('You have input wrong arguments! \n');
    return;
end



