%% User settings

WINLIB1 = [matlabroot '\extern\lib\win32\lcc\libmwlapack.lib'];
WINLIB2 = [matlabroot '\extern\lib\win32\lcc\libmwblas.lib'];

%% Code

% if strcmp(computer,'PCWIN64') | strcmp(computer,'GLNXA64')
%   computer_model = 64; 
% else
%   computer_model = 32; 
% end

matlabversion = sscanf(version,'%f');
matlabversion = matlabversion(1);

% fsp = filesep;

% subdir = strcat('source',fsp);
subdir = 'source/';

fname{ 1} = 'mexsdplr.c';
fname{ 2} = 'copystructures.c';
fname{ 3} = 'dataoper.c';
fname{ 4} = 'eigen.c';
fname{ 5} = 'initialize.c';
fname{ 6} = 'lbfgs.c';
fname{ 7} = 'linesearch.c';
fname{ 8} = 'main.c';
fname{ 9} = 'misc.c';
fname{10} = 'params.c';
fname{11} = 'rankreduce.c';
fname{12} = 'readdata.c';
fname{13} = 'sdplrlib.c';
fname{14} = 'timefuncs.c';
fname{15} = 'util.c';

mexcmd = 'mex -g -O -v CFLAGS="\$CFLAGS -std=iso9899:1999" -D__MEX ';
%mexcmd = [mexcmd '-D__MEXSS '];

if ispc
 mexcmd = [mexcmd '-D__WIN32 '];
end
if matlabversion >= 7.3
 mexcmd = [mexcmd '-largeArrayDims ']; 
end
for k = 1:length(fname)
 mexcmd = [mexcmd subdir fname{k} ' '];
end 
mexcmd = [mexcmd 'gsl-1.5/poly/eval.c gsl-1.5/poly/solve_cubic.c '];
if ispc
 mexcmd = [mexcmd '"' WINLIB1 '"' ' ' '"' WINLIB2 '"' ];
else
 %mexcmd = [mexcmd 'libarpack_LINUX.a -lmwlapack -lmwblas libg2c.so'];
 %%uncomment the above line if using ARPACK
 mexcmd = [mexcmd '-lmwlapack -lmwblas'];
end

mexcmd
eval(mexcmd)
