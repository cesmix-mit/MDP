%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

% Add MDP to Matlab search path
addpath(char(MDPpath + "src/Matlab/Gencode"));
addpath(char(MDPpath + "src/Matlab/Preprocessing"));
addpath(char(MDPpath + "src/Matlab/Postprocessing"));
addpath(char(MDPpath + "src/Matlab/Fitting"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');

