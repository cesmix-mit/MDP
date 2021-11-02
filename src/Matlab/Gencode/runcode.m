%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function runcode(app)

bindir = "exec";
cd(char(app.sourcepath + bindir));
if app.platform == "cpu"
    eval("!./cpuMDP " + app.appname + " out");
elseif app.platform == "gpu"
    eval("!./gpuMDP " + app.appname + " out");    
end
cd(char(app.currentdir));

end

