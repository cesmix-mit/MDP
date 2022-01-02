
__precompile__()

module MDP

include("Preprocessing/Preprocessing.jl");
include("Gencode/Gencode.jl");
include("Postprocessing/Postprocessing.jl");

export initializeapp, initializeconfig, initializemdp, preprocessing, readconfig, readEXTXYZ, readJSON, readJSONfromFolder
export setlattice, setregion, setdomain

end


