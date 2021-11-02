%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [app,config] = initializemdp(sourcepath,version)

app = initializeapp(sourcepath,version);
config = initializeconfig(app);

end
