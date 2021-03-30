function cppfiles(foldername,filename, stropu, strcpu, strgpu, potnum)

#foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

if potnum==1
    stropu = replace(stropu, "int ng)" => "int ng, int potnum)");
    stropu = replace(stropu, "int);" => "int, int);");
    strcpu = replace(strcpu, "int ng)" => "int ng, int potnum)");
    strcpu = replace(strcpu, "int);" => "int, int);");
    strgpu = replace(strgpu, "int ng)" => "int ng, int potnum)");        
    strgpu = replace(strgpu, "int);" => "int, int);");
end

ioopu = open(foldername * "/" * opufile * ".cpp", "w");
write(ioopu, stropu);
close(ioopu);

iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
write(iocpu, strcpu);
close(iocpu);

iogpu = open(foldername * "/" * gpufile * ".cu", "w");
write(iogpu, strgpu);
close(iogpu);

end

