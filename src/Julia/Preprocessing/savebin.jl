function savebin(filename,a)
    fileID = open(filename,"w");
    write(fileID,Float64.(a[:]));
    close(fileID);
end

