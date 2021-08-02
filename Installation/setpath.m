%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

% Add Exasim to Matlab search path
addpath(char(sourcepath + "Matlab/Gencode"));
addpath(char(sourcepath + "Matlab/Preprocessing"));
addpath(char(sourcepath + "Matlab/Postprocessing"));
addpath(char(sourcepath + "C++/Main"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');

