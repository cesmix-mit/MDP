function cppfiles(filename, stropu, strcpu, strgpu, potnum)

foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

if potnum==1
    stropu = strrep(stropu, "int ng)", "int ng, int potnum)");
    stropu = strrep(stropu, "int);", "int, int);");
    strcpu = strrep(strcpu, "int ng)", "int ng, int potnum)");
    strcpu = strrep(strcpu, "int);", "int, int);");
    strgpu = strrep(strgpu, "int ng)", "int ng, int potnum)");        
    strgpu = strrep(strgpu, "int);", "int, int);");
end

fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

fid = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(fid);

fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(fid);    

end

