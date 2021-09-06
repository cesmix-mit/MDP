%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

% Add MDP to Matlab search path
addpath(char(MDPpath + "Matlab/Gencode"));
addpath(char(MDPpath + "Matlab/Preprocessing"));
addpath(char(MDPpath + "Matlab/Postprocessing"));
addpath(char(MDPpath + "C++/Main"));
addpath(char(MDPpath + "Installation"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');

