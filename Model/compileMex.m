% compileMex.m — compile dPUdT_core_mex.c into a MATLAB MEX binary.
%
%   Run this script once (or whenever dPUdT_core_mex.c or packMexVals.m
%   is changed) to produce dPUdT_core_mex.mexw64.
%
%   Uses GCC from OpenModelica's bundled MinGW-w64 (ucrt64) since no
%   standalone MinGW Add-On is installed. Writes a temporary batch file
%   to handle Windows paths-with-spaces correctly.
%
%   After compilation, enable the MEX path by adding to your params script:
%     params0.UseFastPPval   = true;
%     params0.UseCompiledMex = true;
%     params0.modelFcn = '@dPUdT_CombinedTransitions_mex';
%
%   See also: packMexVals, dPUdT_CombinedTransitions_mex

model_dir = fileparts(mfilename('fullpath'));
src     = fullfile(model_dir, 'dPUdT_core_mex.c');
out_mex = fullfile(model_dir, 'dPUdT_core_mex.mexw64');
ml      = matlabroot;
lib_dir = fullfile(ml, 'extern', 'lib', 'win64', 'mingw64');
expdef  = fullfile(lib_dir, 'mexFunction.def');
inc_dir = fullfile(ml, 'extern', 'include');

% GCC from OpenModelica — change this path if using a different MinGW
gcc_bin = 'C:\Program Files\OpenModelica1.26.1-64bit\tools\msys\ucrt64\bin';

bat_file = [tempdir 'compile_dPUdT_mex.bat'];
fid = fopen(bat_file, 'w');
fprintf(fid, '@echo off\r\n');
fprintf(fid, 'set PATH=%s;%%PATH%%\r\n', gcc_bin);
fprintf(fid, 'gcc -O2 -DNDEBUG -DMATLAB_MEX_FILE -m64 -shared ^\r\n');
fprintf(fid, '  -fexceptions -fno-omit-frame-pointer ^\r\n');
fprintf(fid, '  -Wl,--no-undefined ^\r\n');
fprintf(fid, '  "-Wl,%s" ^\r\n', expdef);
fprintf(fid, '  "-I%s" ^\r\n', inc_dir);
fprintf(fid, '  "%s" ^\r\n', src);
fprintf(fid, '  "%s\\libmx.lib" ^\r\n', lib_dir);
fprintf(fid, '  "%s\\libmex.lib" ^\r\n', lib_dir);
fprintf(fid, '  "%s\\libmat.lib" -lm ^\r\n', lib_dir);
fprintf(fid, '  -o "%s"\r\n', out_mex);
fclose(fid);

fprintf('Compiling %s ...\n', src);
[status, result] = system(bat_file);
if ~isempty(result)
    fprintf('%s\n', result);
end
delete(bat_file);

if status == 0 && exist(out_mex, 'file')
    fprintf('Done: %s\n', out_mex);
else
    error('compileMex: GCC failed (exit %d). Check GCC path in compileMex.m.', status);
end
