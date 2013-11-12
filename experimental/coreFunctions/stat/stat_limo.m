function stat_limo(varargin,nboot)

% calls functions from the limo_eeg_toolbox to do stats in SIFT
% Cyril Pernet, Tim Mullen - around August/Septembre 2013

% nboot = number of bootstrap to do

% stat_limo(1,data,nboot);
%                    1                = a one-sample t-test
%                    y                = data (dim node, node, freq, time, subjects)
%                   
% stat_limo(2,data1,data2,nboot);
%                    2                = two samples t-test
%                    y1               = data (dim node, node, freq, time, subjects)
%                    y2               = data (dim node, node, freq, time, subjects)
%
% stat_limo(3,data1,data2,nboot);
%                    3                = paired t-test
%                    y1               = data (dim node, node, freq, time, subjects)
%                    y2               = data (dim node, node, freq, time, subjects)
%
% stat_limo(4,y,X,parameter number,nboot);
%                    4                = regression analysis
%                    y                = data (dim node, node, freq, time, subjects)
%                    X                = continuous regressor(s)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%
% stat_limo(5,y,cat,cont,nboot);
%                    5     = N-way ANOVA/ANCOVA
%                    y     = data (dim node, node, freq, time, subjects)
%                    cat   = categorical variable(s) ie groups eg [12121212];
%                    cont  = continuous regressors (covariates) 
%
% stat_limo(6,y,gp,factor_levels,nboot);
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim node, node, freq, time, subjects, measures)
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor


%% do some input checks

type = varargin{1};

% these could be passed as input as well - name would be condition1 for
% instance
name = 'test';
boot_name = ['boot_' name];


%% run the stats

switch type

    %--------------------------------------------------------------------------
    %            One Sample t-test     //     bootstrap-t method
    %--------------------------------------------------------------------------
    case {1}
        
        [n,m,f,t,s]=size(data);
        % ------------------------------------------------
        % make a one_sample file (n*m*freq*time*[t, p])
        % se and df can be used to build CI for display
        one_sample = NaN(n,m,t,f,2);   
        for edge1 = 1:n
            for edge2 =1:m
            Y = squeeze(data(m,n,:,:,find(~isnan(data(m,n,:,:,:))))); % time*freq minus missing guys
            [one_sample(n,m,:,:,1),tm,tci,se,one_sample(n,m,:,:,2),tc,df]=limo_trimci(Y);
            clear tmp Y
            end
        end
        save ([name],'one_sample', '-v7.3')

        % ------------------------------------------------
        % Bootstrap
        % create a boot one_sample file to store data under H0 
        H0_one_sample = NaN(m,n,t,f,2,nboot); 
        centered_data = data - repmat(nanmean(data,5),[1 1 1 1 s]);
        boot_table = randi(s,s,599);
        
        % get results under H0 
        for edge1 = 1:n
            for edge2 =1:m
                for b=1:nboot
                    tmp = squeeze(data(m,n,:,:,boot_table(:,b)));
                    Y = squeeze(tmp(m,n,:,:,find(~isnan(tmp(m,n,:,:,:))))); % time*freq minus missing guys
                    [H0_one_sample(n,m,:,:,1,b),tm,tci,se,H0_one_sample(n,m,:,:,2,b),tc,df]=limo_trimci(Y);
                    clear tmp Y
                end
            end
        end
        save ([boot_name],'H0_one_sample','-v7.3');
            

        
        
        %--------------------------------------------------------------------------
        %            Two Samples t-test     //     percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {2}
        
        [n,m,f,t,s1]=size(data1);
        [n,m,f,t,s2]=size(data2);
        % ------------------------------------------------ 
        % make a two_samples file per parameter (n*m*freq*time*[t, p])
        two_samples = NaN(n,m,f,t,2);
        
        for edge1 = 1:n
            for edge2 =1:m
                Y1 = squeeze(data1(n,m,:,:,find(~isnan(data1(n,m,:,:,:)))));
                Y2 = squeeze(data2(n,m,:,:,find(~isnan(data2(n,m,:,:,:)))));
                [two_samples(n,m,:,:,1),d,se,ci,two_samples(n,m,:,:,2),tc,df]=limo_yuen_ttest(Y1,Y2); clear Y1 Y2
            end
        end
        save ([name],'two_samples', '-v7.3')
        
        % ------------------------------------------------
        % create a boot one_sample file to store data under H0
        H0_two_samples = NaN(n,n,t,f, 2, nboot); 
        data1_centered = data1 - repmat(nanmean(data1,5),[1 1 1 1 s1]);
        data2_centered = data2 - repmat(nanmean(data2,5),[1 1 1 1 s2]);
        boot_table1 = randi(s1,s1,599);
        boot_table2 = randi(s2,s2,599);
        
        % get results under H0
        for edge1 = 1:n
            for edge2 =1:m
                for b=1:nboot
                    tmp = squeeze(data1(m,n,:,:,boot_table1(:,b)));
                    Y1 = squeeze(tmp(m,n,:,:,find(~isnan(tmp(m,n,:,:,:))))); % time*freq minus missing guys
                    clear tmp; tmp = squeeze(data2(m,n,:,:,boot_table2(:,b)));
                    Y2 = squeeze(tmp(m,n,:,:,find(~isnan(tmp(m,n,:,:,:))))); % time*freq minus missing guys
                    [H0_two_samples(n,m,:,:,1,b),d,se,ci,H0_two_samples(n,m,:,:,2,b),tc,df]=limo_yuen_ttest(Y1,Y2); clear Y1 Y2
                end
            end
        end
        save ([boot_name],'H0_two_samples','-v7.3');
            
        
        
        %--------------------------------------------------------------------------
        %            Paired t-test     //     percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {3}
        
        [n,m,f,t,s1]=size(data1);
        % ------------------------------------------------
        % make a paired_samples file per parameter (n*m*freq*time*[t, p])
        paired_samples = NaN(n,m,t,p,2);
        for edge1 = 1:n
            for edge2 =1:m
                Y1 = squeeze(data1(n,m,:,:,find(~isnan(data1(n,m,:,:,:)))); 
                Y2 = squeeze(data1(n,m,:,:,find(~isnan(data1(n,m,:,:,:)))); 
                [paired_samples(n,m,:,:,1),diff,se,ci,paired_samples(n,m,:,:,2),tc,df]=limo_yuend_ttest(Y1,Y2); clear Y1 Y2
            end
        end
        save ([name],'paired_samples', '-v7.3')
        
        % ------------------------------------------------
        % create a boot one_sample file to store data under H0
        H0_paired_samples = NaN(n,m,t,f,2, nboot); % stores T and p values for each boot
        % create centered data to estimate H0
        data1_centered = data1 - repmat(nanmean(data1,5),[1 1 1 1 s]);
        data2_centered = data2 - repmat(nanmean(data2,5),[1 1 1 1 s]);
        boot_table = randi(s,s,599);
        
        % get results under H0
        for edge1 = 1:n
            for edge2 =1:m
                for b=1:nboot
                    tmp = squeeze(data1(m,n,:,:,boot_table(:,b)));
                    Y1 = squeeze(tmp(m,n,:,:,find(~isnan(tmp(m,n,:,:,:))))); % time*freq minus missing guys
                    clear tmp; tmp = squeeze(data2(m,n,:,:,boot_table(:,b)));
                    Y2 = squeeze(tmp(m,n,:,:,find(~isnan(tmp(m,n,:,:,:))))); % time*freq minus missing guys
                    [H0_paired_samples(n,m,:,:,1,b),diff,se,ci,H0_paired_samples(n,m,:,:,2,b),tc,df]=limo_yuend_ttest(Y1,Y2); clear Y1 Y2
                end
            end
        end
        save ([boot_name],'H0_paired_samples','-v7.3');
        
        
        %--------------------------------------------------------------------------
        %                    Regression  
        %--------------------------------------------------------------------------
    case {4}
        
        [n,m,f,t,s]=size(data);
        if s~=size(X,1) % matrix of covariates
            error('the regressor matrix is of the wrong dimension');
        end
        
        % ------------------------------------------------
        % make Regressors file to save the data dim(n*m*freq*time*number of regressors*[F p])
        GLM(n,m,f,t,size(X,2),2)
        H0_GLM(n,m,f,t,size(X,2),2,nboot)
        Boot_table = randi(s,s,599);
        X = [zscore(X) ones(1,sized(X,1))];
        
        % the time*freq being not to big we could make a vector f*t, pass this and reshape avoiding a loop
        % but shitdim is slow as well ; need testing what is faster
        % data = shitfdim(reshape(shitfdim(data,4)),[s,n,m,f*t]),1);
        
        for edge1 = 1:n
            for edge2 =1:m
                for freq =1:f  
                    model = limo_glm1(squeeze(data(n,m,f,:,:))',X,0,0,size(regressor,2),'OLS');
                    GLM(n,m,:,f,1,:) = model.conditions.F;
                    GLM(n,m,:,f,2,:) = model.conditions.p;
                    
                    % boot data
                    model = limo_glm1_boot(squeeze(data(n,m,f,:,:))',X,0,0,size(regressor,2),Z,'OLS',Boot_table);
                    H0_GLM(n,m,:,f,1,:,:) = cell2mat(model.conditions.F);
                    H0_GLM(n,m,:,f,2,:,:) = cell2mat(model.conditions.p);
                end
            end
        end
        save ([name],'GLM','-v7.3');
        save ([boot_name],'H0_GLM','-v7.3');

        
        %--------------------------------------------------------------------------
        %                    N-ways ANOVA / ANCOVA
        %--------------------------------------------------------------------------
    case {5}
        
        [n,m,f,t,s]=size(data);
        if s~=size(Cat,1) % matrix of groups
            error('the condition matrix is of the wrong dimension');
        end
        if ~isempty(Cont)
            if s~=size(Cont,1) % matrix of covariates
            error('the covariate matrix is of the wrong dimension');
            end
        end
        
        % because we need a design matrix - of dim s*nb_conditions it's
        % easier to reshape the data apply the GLM and put them back into
        % original format - but shitdim is slow
        % ------------------------------------------------
        Y = shiftdim(reshape(shiftdim(data,2),[f,t,s,n*m]),3);
        clear data; Boot_table = randi(s,s,599);

        % make design matrix and files
        if isempty(Cont)
            [X,nb_conditions,nb_interactions,nb_continuous] = limo_design_matrix(squeeze(Y(:,1,:,:)),Cat,Cont,pwd,1,0,1);
        else
            [X,nb_conditions,nb_interactions,nb_continuous] = limo_design_matrix_tf(Y,Cat,Cont,pwd,1,1,1); 
        end
        
        delete Betas.mat R2.mat Res.mat Yhat.mat Yr.mat % stuff created by limo_design_matrix
        a = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
         if strcmp(a,'No')
             return
         else
             GLM = NaN(n*m,f,t,size(Cat,2),2); % store F, p  values 
             H0_GLM = NaN(n*m,f,t,size(Cat,2),2,nboot); % store F, p  values 
             for edge = 1:(n*m)
                 fprintf('computing anova edge %g/%g',edge,n*m); disp(' ');
                     for freq =1:f  %% the time*freq being not to big we could make then vector, pass this and reshape avoiding this loop
                         % observed data
                         model = limo_glm1(squeeze(Y(edge,f,:,:))',X,nb_conditions,nb_interactions,nb_continuous,'OLS');
                         GLM(edge,f,:,:,1) = model.conditions.F;
                         GLM(edge,f,:,:,2) = model.conditions.p;
                         
                         % boot data
                         model = limo_glm1_boot(squeeze(Y(edge,f,:,:))',X,nb_conditions,nb_interactions,nb_continuous,1,'OLS',Boot_table);
                         H0_GLM(edge,f,:,:,1,:) = cell2mat(model.conditions.F);
                         H0_GLM(edge,f,:,:,2,:) = cell2mat(model.conditions.p);
                     end
             end
              % need to reshape to standard dim
              save ([name],'GLM','-v7.3');
              save ([boot_name],'H0_GLM','-v7.3');
         end
         
        
        
       %----------------------------------------------------------------------------------------------
        %                    Repeated Measure ANOVA (multivariate approach) - bootstrap centering data
        %----------------------------------------------------------------------------------------------
    case {6}
        
        [n,m,f,t,s]=size(data);
        
        % specific stuff for repeated measures
        % from the input we know which case to handle
        if unique(gp_vector) == 1
            % one sample
            if length(factor_levels) ==1
                type = 1;
            elseif length(factor_levels) >1
                type = 2;
            end
        else
            % k samples
            if length(factor_levels) ==1
                type = 3;
            elseif length(factor_levels) >1
                type = 4;
            end
        end
        
         
        % make files to be stored
        if type == 1 % one factor
            C = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            X = [kron(eye(prod(factor_levels)),ones(size(data,3),1)) ones(size(x,1),1)];
            tmp_Rep_ANOVA = NaN(m,n,f,t,1,2); % store F and p (means can computed via the central tendency and CI button)
            
        elseif type == 2 % many factor
            C = limo_OrthogContrasts(factor_levels);
            index = length(factor_levels)+1;
            X = [kron(eye(prod(factor_levels)),ones(size(data,3),1)) ones(size(x,1),1)];
            tmp_Rep_ANOVA = NaN(m,n,f,t,length(C),2); % store F and p for each within factor and interactions
            
        elseif type == 3 % one factor within and one factor between
            gp_values = unique(gp_vector); 
            k = length(gp_values); 
            X = NaN(size(gp_vector,1),k+1);
            for g =1:k 
                X(:,g) = gp_vector == gp_values(g); 
            end 
            X(:,end) = 1; % design matrix for gp effects
            C  = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            X = [kron(X(:,1:k),eye(prod(factor_levels))) sum(x,2)]; % just for display
            tmp_Rep_ANOVA                     = NaN(m,n,f,t,1,2);
            Rep_ANOVA_Gp_effect               = NaN(m,n,f,t,2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(m,n,f,t,1,2);
             
        elseif type == 4 % many factors within and one factor between
            gp_values = unique(gp_vector); 
            k = length(gp_values); 
            X = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g); 
            end
            X(:,end) = 1; % design matrix for gp effects
            C = limo_OrthogContrasts(factor_levels);
            X = [kron(X(:,1:k),eye(prod(factor_levels))) sum(x,2)]; % just for display
            tmp_Rep_ANOVA                     = NaN(m,n,f,t,length(C),2);
            Rep_ANOVA_Gp_effect               = NaN(m,n,f,t,2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(m,n,f,t,length(C),2);
       end
        
        % check the design with user
        figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc(X);
        colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors');
        ylabel('subjects'); drawnow;
        go = questdlg('start the analysis?');
        if strcmp(go,'No')
            return
        end
        
        % do the analysis
        index = 1;
        for edge1 = 1:n
            for edge2 =1:m
            fprintf('analyse edge %g/%g\n ...', edge,n*m); 
            index= index+1;
            
            tmp = squeeze(data(edge1,edge2,:,:,:,:));
            Y = tmp(:,:,find(~isnan(tmp(1,1,:,1))),:);
            gp = gp_vector(find(~isnan(tmp(1,1,:,1))),:);
            if type == 3 || type == 4
                XB = X(find(~isnan(tmp(1,1,:,1))),:);
            end
            
            if type == 1
                for freq = 1:f
                    result = limo_rep_anova(squeeze(Y(f,:,:,:)),gp,factor_levels,C); % usual means
                    tmp_Rep_ANOVA(edge1,edge2,f,:,1,1) = result.F;
                    tmp_Rep_ANOVA(edge1,edge2,f,:,1,2) = result.p;
                end
                
            elseif type == 2
                for freq = 1:f
                    result = limo_rep_anova(squeeze(Y(f,:,:,:)),gp,factor_levels,C); % usual means
                    tmp_Rep_ANOVA(edge1,edge2,f,:,:,1) = result.F';
                    tmp_Rep_ANOVA(edge1,edge2,f,:,:,2) = result.p';
                end
                
            elseif type == 3
                for freq = 1:f
                    result = limo_rep_anova(squeeze(Y(f,:,:,:)),gp,factor_levels,C,XB); % usual means
                    tmp_Rep_ANOVA(edge1,edge2,f,:,1,1)                     = result.repeated_measure.F;
                    tmp_Rep_ANOVA(edge1,edge2,f,:,1,2)                     = result.repeated_measure.p;
                    Rep_ANOVA_Gp_effect(edge1,edge2,f,:,1)                 = result.gp.F;
                    Rep_ANOVA_Gp_effect(edge1,edge2,f,:,2)                 = result.gp.p;
                    tmp_Rep_ANOVA_Interaction_with_gp(edge1,edge2,f,:,1)   = result.interaction.F;
                    tmp_Rep_ANOVA_Interaction_with_gp(edge1,edge2,f,:,2)   = result.interaction.p;
                end
                
            elseif type == 4
                for freq = 1:f
                    result = limo_rep_anova(squeeze(Y(f,:,:,:)),gp,factor_levels,C,XB); % usual means
                    tmp_Rep_ANOVA(edge1,edge2,f,:,:,1)                     = result.repeated_measure.F';
                    tmp_Rep_ANOVA(edge1,edge2,f,:,:,2)                     = result.repeated_measure.p';
                    Rep_ANOVA_Gp_effect(edge1,edge2,f,:,1)                 = result.gp.F;
                    Rep_ANOVA_Gp_effect(edge1,edge2,f,:,2)                 = result.gp.p;
                    tmp_Rep_ANOVA_Interaction_with_gp(edge1,edge2,f,:,:,1) = result.interaction.F';
                    tmp_Rep_ANOVA_Interaction_with_gp(edge1,edge2,f,:,:,2) = result.interaction.p';
                end
            end
            clear tmp Y gp result
        end
        
        
        % save stuff
        % ---------
        for i=1:size(tmp_Rep_ANOVA,3)
            name = sprintf('Rep_ANOVA_Factor_%g',i);
            % save each factor effect as F/p values
            Rep_ANOVA = squeeze(tmp_Rep_ANOVA(:,:,:,:,i,:)); 
            save([name],'Rep_ANOVA', '-v7.3'); 
        end
        
        if type == 3 || type ==4
            for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
                name = sprintf('Rep_ANOVA_Interaction_gp_Factor_%g',i);
                % save each interaction effect as F/p values
                Rep_ANOVA_Interaction_with_gp = squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,:,:,i,:));
                save([name],'Rep_ANOVA_Interaction_with_gp', '-v7.3'); clear Rep_ANOVA_Interaction_with_gp;
            end
            save Rep_ANOVA_Gp_effect Rep_ANOVA_Gp_effect -v7.3; % always only 1 effect
        end
        
        
        % ----------------------------------------------------------------
%         if nboot > 0
%             
%                 % create files to store bootstrap under H1 and H0
%                 disp('making bootstrap files ...')
%                 if type ==1
%                     tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2,nboot);
%                 elseif type == 2
%                     tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2,nboot);
%                 elseif type == 3
%                     tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2,nboot);
%                     boot_H0_Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2,nboot);
%                     tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),1,2,nboot);
%                 else
%                     tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2,nboot);
%                     boot_H0_Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2,nboot);
%                     tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),length(C),2,nboot);
%                 end
%             
%                 % the data have to be centered (H0) for each cell
%                 centered_data = NaN(size(data,1),size(data,2),size(data,3),size(data,4));           
%                 nb_conditions = prod(factor_levels);
% 
%                 if type ==1 || type ==2
%                     for condition=1:nb_conditions
%                         avg = repmat(nanmean(data(:,:,:,condition),3),[1 1 size(data,3)]);
%                         centered_data(:,:,:,condition) = data(:,:,:,condition) - avg; clear avg
%                     end
%                 else
%                     for gp=1:LIMO.design.nb_conditions % = nb gp
%                         gp_index = find(LIMO.data.Cat == gp);
%                         for condition=1:nb_conditions
%                             avg = repmat(nanmean(data(:,:,gp_index,condition),3),[1 1 length(gp_index)]);
%                             centered_data(:,:,gp_index,condition) = data(:,:,gp_index,condition) - avg; clear avg
%                         end
%                     end
%                 end
%                 save(['H0', filesep, 'centered_data'], 'centered_data', '-v7.3');
% 
%                 % create an index to use across all electrodes and frames
%                 % (different per gp but identical across conditions)
%                 disp('making random table...')
%                 % boot_table = limo_create_boot_table(squeeze(data(:,:,:,1)),nboot);
%                 boot_table = limo_create_boot_table(reshape(data(:,:,:,1),[size(data,1) size(data,2) size(data,3)]),nboot);
%                 save(['H0', filesep, 'boot_table'], 'boot_table')
%             
%                 % compute bootstrap under H0 for F and p
%                 for B=1:nboot
%                     fprintf('Repeated Measures ANOVA bootstrap %g \n ...', B);
%                     for e = 1:length(array)
%                         electrode = array(e);
%                         
%                         tmp = squeeze(centered_data(electrode,:,boot_table{electrode}(:,B),:));
%                         Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
%                         gp = gp_vector(find(~isnan(tmp(1,:,1))));
%                         if type == 3 || type == 4
%                             XB = X(find(~isnan(tmp(1,:,1))));
%                         end
%                         
%                         if type == 1
%                             if strcmp(LIMO.design.method,'TM')
%                                 result = limo_robust_rep_anova(Y,gp,factor_levels,C);
%                             else
%                                 result = limo_rep_anova(Y,gp,factor_levels,C);
%                             end
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,1,1,B)         = result.F;
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,1,2,B)         = result.p;
%                         elseif type == 2
%                             if strcmp(LIMO.design.method,'TM')
%                                 result = limo_robust_rep_anova(Y,gp,factor_levels,C);
%                             else
%                                 result = limo_rep_anova(Y,gp,factor_levels,C);
%                             end
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,:,1,B)         = result.F';
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,:,2,B)         = result.p';
%                         elseif type == 3
%                             if strcmp(LIMO.design.method,'TM')
%                                 result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
%                             else
%                                 result = limo_rep_anova(Y,gp,factor_levels,C,XB);
%                             end
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,1,1,B)                     = result.repeated_measure.F;
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,1,2,B)                     = result.repeated_measure.p;
%                             H0_Rep_ANOVA_Gp_effect(electrode,:,1,B)                      = result.gp.F;
%                             H0_Rep_ANOVA_Gp_effect(electrode,:,2,B)                      = result.gp.p;
%                             tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,1,1,B) = result.interaction.F;
%                             tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,1,2,B) = result.interaction.p;
%                         elseif type == 4
%                             if strcmp(LIMO.design.method,'TM')
%                                 result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
%                             else
%                                 result = limo_rep_anova(Y,gp,factor_levels,C,XB);
%                             end
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,:,1,B)                     = result.repeated_measure.F';
%                             tmp_boot_H0_Rep_ANOVA(electrode,:,:,2,B)                     = result.repeated_measure.p';
%                             H0_Rep_ANOVA_Gp_effect(electrode,:,1,B)                      = result.gp.F;
%                             H0_Rep_ANOVA_Gp_effect(electrode,:,2,B)                      = result.gp.p;
%                             tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,:,1,B) = result.interaction.F';
%                             tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,:,2,B) = result.interaction.p';
%                             clear y result
%                         end
%                         clear XB Y gp tmp
%                     end
%                 end
%                 
%                 % save 
%                 for i=1:size(tmp_boot_H0_Rep_ANOVA,3)
%                     name = sprintf('H0_Rep_ANOVA_Factor_%g',i);
%                     H0_Rep_ANOVA = squeeze(tmp_boot_H0_Rep_ANOVA(:,:,i,:,:)); % save each factor effect as F/p/nboot values
%                     save(['H0', filesep, name],'H0_Rep_ANOVA', '-v7.3'); 
%                 end
%                 
%                 if type == 3 || type ==4
%                     for i=1:size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp,3)
%                         name = sprintf('H0_Rep_ANOVA_Interaction_gp_Factor_%g',i);
%                         H0_Rep_ANOVA_Interaction_with_gp = squeeze(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,i,:)); % save each interaction effect as F/p values
%                         save(['H0', filesep, name],'H0_Rep_ANOVA_Interaction_with_gp', '-v7.3'); clear H0_Rep_ANOVA_Interaction_with_gp;
%                     end
%                     
%                     save(['H0', filesep, 'H0_Rep_ANOVA_Gp_effect'], 'H0_Rep_ANOVA_Gp_effect', '-v7.3'); 
%                 end
%             end
            
end


