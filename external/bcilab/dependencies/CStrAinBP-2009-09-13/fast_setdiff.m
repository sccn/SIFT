function A = fast_setdiff(A,B)
% A fast version of setdiff for cell arrays of strings.

try
    A(CStrAinBP(A,B)) = [];
catch
    % ... and the code to compile it, if necessary
    persistent action_taken; %#ok<TLEV>
    % fall back to regular setdiff
    A = setdiff(A,B);

    % try to fix the problem
    if isempty(action_taken)
        e = lasterror; %#ok<LERR>
        try
            % try to recompile it ad hoc
            dir = pwd;
            c = onCleanup(@()cd(dir));
            basepath = fileparts(mfilename('fullpath'));
            cd(basepath);
            if ~exist('noretry.txt','file')
                disp('BCILAB can use a 20x faster "setdiff" function if it succeeds to compile it for your platform. Trying now...');
                if ispc && strfind(computer,'64')
                    disp('Note that, since you are using 64-bit Windows, this can be challenging.');
                    disp('You probably need to install the (free) Microsoft Visual Studio 2005 Express (or optionally 2008 Express for MATLAB 2008b+');
                    disp('or 2010 Express for MATLAB 2010a+) and the Microsoft Platform SDK (6.1 for 2008, 7.1 for 2010) for your Windows Version.');
                    disp('See also: http://argus-home.coas.oregonstate.edu/forums/development/core-software-development/compiling-64-bit-mex-files');
                    disp('          http://www.mathworks.com/support/compilers/R2010b/win64.html');
                elseif isunix
                    disp('Since you are using Linux/UNIX, you are probably using your system''s GCC compiler by default.');
                    disp('If the installation fails, consider installing a version of GCC that is compatible with your MATLAB version,');
                    disp('or fall back to the builtin LCC compiler that ships with MATLAB.');
                    disp('Supported GCC versions include: R14SP3->3.2(x64)/3.3(x32), R2006a->3.4.4, R2006b->3.4.5, R2007a->4.1.1');
                    disp('R2007b->4.0-4.2, R2008a->4.0-4.2, R2008b->4.0-4.2, R2009a->4.1-4.2, R2009b->4.2.3, R2010a->4.2.3, R2010b->4.3.x');
                end
                try
                    mex CStrAinBP.cpp % try the C++ version (gcc, etc.)
                catch
                    mex CStrAinBP.c   % fall back to the C version (lcc, etc.)
                end
                if isequal(2,CStrAinBP({'a','b'},{'b','c'}))
                    % successful: put it in a proper subfolder
                    [dum,host] = system('hostname'); host = strtrim(host);%#ok<ASGLU>
                    io_mkdirs([basepath filesep host filesep],{'+w','a'});
                    if ispc
                        system(['move CStrAinBP.mex* ' host filesep]);
                        system(['copy env_add.reference ' host filesep 'env_add.m']);
                    else
                        system(['mv CStrAinBP.mex* .' filesep host filesep]);
                        system(['cp env_add.reference .' filesep host filesep 'env_add.m']);
                    end
                    addpath([basepath filesep host filesep]);
                    disp('Successfully recompiled fast_setdiff()''s backend function for your machine.');
                else
                    % failed: delete the file...
                    if ispc
                        system('del CStrAinBP.mex*');
                    else
                        system('rm CStrAinBP.mex*');
                    end
                    error('Compilation failed.');
                end
            end
        catch
            warning('BCILAB:slow_setdiff','Falling back to the regular setdiff() implementation, as no compatible binary could be created at this point.');
            f = fopen('noretry.txt','w+');
            fprintf(f,'Delete this file if you want to retry auto-compilation after you have chosen a compiler via mex -setup');
            fclose(f);
        end
        action_taken = true;
    end
end