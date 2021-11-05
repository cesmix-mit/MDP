cdir = pwd(); ii = strfind(cdir, "MDP");
if isempty(ii) == 0
    MDPpath = cdir(1:(ii+2)) + "/";
    run(MDPpath + "installation/setpath.m");
else
    % MDPpath = /path/to/MDP;
    disp("MDP's path is not found. Please uncomment the above line and set the path to MDP.");
end

species = ["Ta"];
datapath = "/JSON/Displaced_A15";
config = readJSONfromFolder(datapath, species);

