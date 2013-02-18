function [S,W,E] = FastICA(Z,varargin)
% FastICA ---- Fast Fixed-Point ICA algorithm
%     由输入的白化后信号（零均值单位方差就）以及相关参数，采用symmetric正交化或者
%     deflation正交化计算出恢复的原始信号S以及对应的加权矩阵W。当采用symmetric方法时，
%     输入多少个观测信号则估计出多少个独立分量；当采用deflation方法时，可以估计任意数目的独立分量
%     （输出的独立分量数目可小于观测信号数目）
%  
% 函数:
%     [S,W] = FastICA(Z,'method','defl',  'origW',W0,  'NonlinFun','tanh',  'a1',1,  'a2',1, 
%             'maxNumIteration',300,  'OverValue',0.001,  'report','off',  'mixMatrix',VA,  'ICNum',3)
%
% 参数说明......
%       参数名称    ||   参数取值       || 参数说明                   
%               Z                       : 白化了的信号，每行代表一个信号
%         'method'    ['defl'|'symm']   : 采用的正交方法。
%                                         'defl'为deflationary方法
%                                         'symm'为symmetric方法，默认值为'symm'                                                mm'
%          'origW'                      : 输入初始的加权矩阵，
%                                         若省略，则由系统产生初始的正交矩阵
%      'NonlinFunc'    ['pow3'|         : 采用的非线性函数。'pow3'为g(y)=y^3
%                       'tanh'|           'tanh'为g(y)=tanh(a1*y)  
%                       'gaus']           'gaus'为g(y)=y*exp(-y^2/2)                                                            
%              'a1'                     : 非线性函数 'tanh'中调用的值,1<=a1<=2，默认值1
%              'a2'                     : 非线性函数 'gaus'中调用的值，默认值1
% 'maxNumIteration'                     : 最大的迭代次数，默认值为500
%       'OverValue'                     : 决定迭代结束的数值，其越小，表明W*W'越趋于单位阵
%                                         默认值为0.0001
%          'report'    ['on'|'off']     : 是否需要打印出程序运行信息
%                                         'on'为打印，'off'为不打印。默认为'on'
%       'mixMatrix'                     : 即 V*A,V为白化矩阵，A为混合矩阵
%           'ICNum'                     : 输出的估计信号个数，注意此参数仅当'method'参数为'defl'时使用
%                E                      : 记录每次迭代后的性能指数，性能指数由函数CTPI(W'*V*A,'norm')计算
%                S                      : 恢复的信号: S=W'* Z
%                W                      : 恢复信号的加权矩阵,W=[b1,b2,...bm]
%
% 调用格式：
%     如 [S,W]=FastICA(WhitenedSig,'nonlinfunc','POW3','maxNumIteration',100) 
%     输入的可选参数必须成对。先是参数名，接着是对应的参数值。注意输入是数值还是字符串。
%
% 作者：张智林（Zhang Zhi-Lin）
%       zlzhang@uestc.edu.cn          zzl.private@eyou.com
%       http://teacher.uestc.edu.cn/teacher/teacher.jsp?TID=zzl80320
% 
% 版本：1.0    日期：2003年10月31日
% 版本：1.1    更新：2004年12月20日
% 版本：1.2    更新：2005年 7月 1日
% 版本：1.3    更新：2006年 4月11日   （当采用deflation方法时，可估计的信号数目不必等于源信号数目）
%
% 参考文献：
% [1] A.Hyvarinen. Fast and Robust Fixed-Point Algorithms for Independent
%     Component Analysis. IEEE Trans. on Neural Networks,10(3):626-634,1999
%




% ===================================================================
% 取默认值
ICNum=size(Z,1);    % 输入为多少个信号，则恢复出多少个信号    
method='symm';
NonlinFunc='tanh';
a1=1;
a2=1;
maxNumIteration=500;
mixMatrix = eye(ICNum);
OverValue=0.0001;
report='on';
numSamples=size(Z,2);    % 每个信号的采样点数，即信号的长度

UserInputW=0;   % 判断是否由用户输入初始的加权矩阵。为1，则由用户输入。

% ===================================================================
% 获取可选参数
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'method'
                method=varargin{i+1}; method=lower(method);
            case 'origw'  % origW
                W=varargin{i+1}; 
                UserInputW=1;   % 由用户输入初始的加权矩阵
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
                % 输入参数名称有错误
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

if UserInputW == 0
    % 由系统产生初始的正交加权矩阵W
    %W=orth(rand(ICNum)-0.5);
    W=orth(eye(ICNum));
end

% ===================================================================
% 检测输入参数的正确性
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
% 打印各种信息
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
% symmetric 正交方法
if strcmp(method,'symm')
    if( strcmp(report,'on') )   fprintf('Using symmetric method \n'); end
    
    WOld=zeros(size(W));  % WOld用来存贮上一次的W的值，以便和现在的W值相比较
    
    %--------------------------------------------
    % fast fixed-point ICA 算法的迭代
    for loop=1:maxNumIteration + 1
        
        % 若迭代次数到了maxNumIteration次仍没有收敛，则打印信息并返回当前的W，S值。
        if loop == maxNumIteration + 1  
            fprintf('Cannot convergence even after %d steps \n',maxNumIteration);
            S=W'*Z;
            return;
        end
        
        % Symmetric 正交化
        W=W*real(inv(W'*W)^(1/2)); 
        
        % 记录本次迭代的性能指数
        E(loop) = CTPI(W'* mixMatrix,'norm');
        
        % 判断是否符合递归结束的条件。在这里，也考虑到了方向正好相反的向量
        minAbsCos = min(abs(diag(W'*WOld)));  
        
        if(strcmp(report,'on'))
            % 打印迭代过程等信息
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
            return;                %跳出整个循环－－－－－－－－
        end
        
        WOld=W;
        
        % 采用不同的非线性函数计算
        switch NonlinFunc
            % pow3
            case 'pow3'
                W=(Z*((Z'*W).^3))/numSamples-3*W;   % 在整个观测序列上求平均
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
        
        
    end   % 迭代结束
end   % symmetric method 结束


% ===================================================================
% deflation 正交方法,加权矩阵B的每个列向量一个接一个的求出来
if strcmp(method,'defl')
    if( strcmp(report,'on') )   fprintf('Using deflation method \n'); end
    
    B = zeros(size(Z,1),ICNum);   % 最后要求出的加权矩阵
    
    % 处理初始的加权向量
    if UserInputW==0  % 由系统产生初始的加权向量
        w=rand(size(Z,1),1)-0.5;
    else    % 由用户输入初始的加权矩阵
        w=W(:,1);
    end
    w=w-B*B'*w;     % 正交化
    w=w/norm(w);    % 归一化
    wOld=zeros(size(w));   % wOld用来存贮上一次的w值
    
    %--------------------------------------------
    for i=1:ICNum                % 依次求出ICNum个独立分量
        
        % ------------------------------
        % fast fixed-point ICA 算法的迭代
        for loop=1:maxNumIteration+1
            % 已经找到的基向量生成了一个向量空间。现在把当前的向量投影在这个向量空间的
            % 正交空间里。从而实现正交化。注意到之所以能够利用矩阵B正确进行投影，是因
            % 为矩阵B中的零向量对投影不产生任何影响。
            w=w-B*B'*w;     % 正交化
            w=w/norm(w);    % 归一化
            
            if loop == maxNumIteration+1     % 当迭代超过maxNumIteration时，则退出返回
                fprintf('IC No.%d cannot converge in %d iteration.\n',i,maxNumIteration);
                S=B'*Z;
                return;
            end
            
            % 判断是否收敛。当收敛时，w和wOld方向相同，或者正好相反。
            % 所以，判断条件考虑了这两种情况：
            if norm(w-wOld) < OverValue | norm(w+wOld) < OverValue
                if( strcmp(report,'on') )  
                    fprintf('IC No.%d  ---------- Computed after %d iteration\n',i,loop-1);
                end
                B(:,i)=w;  % 将当前的w值保存在B中
                break;
            end
            
            wOld=w;
            
            % 采用不同的非线性函数计算
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
         
        end  % 对B的某个列向量w循环maxNumIteration 次 
       
    end  % 求ICNum个独立分量结束
    
    W=B;     % 加权矩阵的输出变量是W
    S=W'*Z;  % 恢复的信号
    
    if ICNum == size(Z,1)
        E = CTPI(W'* mixMatrix,'norm');   % 性能指数
    end
    
end % deflation method 结束



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



