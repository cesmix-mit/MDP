function genpotnum(filename, npm, sp0, sp1, sp2, sp3)

opufile = "opu" * filename;

tmp = "template <typename T> void " * opufile;
tmp = tmp * sp0;
tmp = tmp * "{\n";
for k = 1:npm
    if k == 1
        tmp = tmp * "\tif (potnum == " * string(k) * ")\n";    
    else            
        tmp = tmp * "\telse if (potnum == " * string(k) * ")\n";    
    end 
    tmp = tmp * "\t\t" * opufile * string(k) * sp3;
end
tmp = tmp * "}\n";
tmp = tmp * "template void " * opufile;
tmp = tmp * sp1;
tmp = tmp * "template void " * opufile;
tmp = tmp * sp2;

stropu = tmp;
strcpu = replace(stropu, "opu" => "cpu");
strgpu = replace(stropu, "opu" => "gpu");

return stropu, strcpu, strgpu

end

