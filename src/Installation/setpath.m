% Add Exasim to Matlab search path
ii = strfind(cdir, "Applications");
versiondir = cdir(1:(ii-1)) + "Matlab";
addpath(char(versiondir + "/Gencode"));
addpath(char(versiondir + "/Preprocessing"));
addpath(char(cdir(1:(ii-1)) + "C++/Main"));

% Set Matlab's PATH enviroment variable so that Exasim can call external packages    
setenv('PATH', '/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin');

