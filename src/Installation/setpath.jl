# Add Exasim to Julia search path
push!(LOAD_PATH, sourcepath * "Julia/Gencode");
push!(LOAD_PATH, sourcepath * "Julia/Preprocessing");
push!(LOAD_PATH, sourcepath * "Julia/Postprocessing");
push!(LOAD_PATH, sourcepath * "C++/Main");

# Set Julia's PATH enviroment variable so that Exasim can call external programs
ENV["PATH"] = "/usr/local/bin:/usr/bin:/opt/local/bin:/bin:/usr/sbin:/sbin";

