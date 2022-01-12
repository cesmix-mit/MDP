function readbin(filename)

    a = reinterpret(Float64,read(filename));
    return a
    
end

