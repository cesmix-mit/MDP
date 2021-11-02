function mstr = gencorestr(mstr)

mstr = strrep(mstr, "template <typename T> T", "inline dstype" );
mstr = strrep(mstr, "template <typename T> void", "inline void" );
mstr = strrep(mstr, "gpu", "" );
mstr = strrep(mstr, "T ", "dstype " );
str0 = mstr;
str0 = strrep(str0, "inline void ", "" );
str0 = strrep(str0, "inline dstype ", "" );
str0 = strrep(str0, "dstype *", "" );
str0 = strrep(str0, "dstype ", "" );
str0 = strrep(str0, "int *", "" );
str0 = strrep(str0, "int ", "" );
mstr = strrep(mstr, ");", ", int backend)" );
mstr = mstr + "\n" + "{" + "\n";
mstr = mstr + "\t" + "if (backend == 1)" + "\n";
mstr = mstr + "\t\t" + "cpu" + str0 + "\n";
mstr = mstr + "#ifdef USE_OMP" + "\n";
mstr = mstr + "\t" + "if (backend == 4)" + "\n";
mstr = mstr + "\t\t" + "omp" + str0 + "\n";
mstr = mstr + "#endif" + "\n";
mstr = mstr + "#ifdef USE_HIP" + "\n";
mstr = mstr + "\t" + "if (backend == 3)" + "\n";
mstr = mstr + "\t\t" + "hip" + str0 + "\n";
mstr = mstr + "#endif" + "\n";
mstr = mstr + "#ifdef USE_CUDA" + "\n";
mstr = mstr + "\t" + "if (backend == 2)" + "\n";
mstr = mstr + "\t\t" + "gpu" + str0 + "\n";
mstr = mstr + "#endif" + "\n";
mstr = mstr + "}" + "\n";

fid = fopen("tmp.cpp", 'w');
fprintf(fid, char(mstr));
fclose(fid);


