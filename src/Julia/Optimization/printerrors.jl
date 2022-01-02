function lengthstring(folders)
    n = 0;
    for i = 1:length(folders)    
        n = max(length(folders[i]), n)    
    end
    return n
end

function replacestring(str1, str2)
     str3 = collect(str1)
     str4 = collect(str2)
     for i = 1:length(str4)
        str3[i] = str4[i]
     end
     return join(str3)
end

function num2string(a, n)
    s = collect(string(a))    
    if length(s) < n
        s2 = repeat('0',(n-length(s)))
        s = join(s) * s2 
    else
        for i = length(s):-1:(n+1) 
            deleteat!(s, i)                           
        end    
        s = join(s)
    end
    return s
end

function printerrors(folders, errors, errstr)

    n = max(lengthstring(folders), length(errstr))
    str1 = repeat(' ', n) * " |   ";

    num = 16;
    print("\n" * repeat('-', 88) * "\n")
    print(replacestring(str1 * "       MAE          |          RMSE          |          RSQ          |\n", errstr))
    print(repeat('-', 88) * "\n")
    e1 = errors[1,1]; e2 = errors[1,2]; e3 = errors[1,3]; 
    print(replacestring(str1 * num2string(e1, num) * "    |    " * num2string(e2, num) * "    |    " * num2string(e3, num) * "  " * " |\n", "ALL"))
    for i = 1:length(folders)
        a1 = errors[i+1,1]; a2 = errors[i+1,2]; a3 = errors[i+1,3]; 
        print(replacestring(str1 * num2string(a1, num) * "    |    " * num2string(a2, num) * "    |    " * num2string(a3, num) * "  " * " |\n", folders[i]))
    end
    print(repeat('-', 88) * "\n")

end
