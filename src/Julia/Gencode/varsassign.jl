function varsassign(mystr::String, varname::String, n::Int, ustr, flg::Int)

if flg==0
    for i = 1:n
        str1 = varname * string(i);
        if contains(ustr, str1)
            str2 = varname * "[" * string(i-1) * "]";
            mystr = mystr * "\t\tT " * str1 * " = " * str2 * ";\n";
        end
    end
elseif flg==1
    for i = 1:n
        str1 = varname * string(i);
        if contains(ustr, str1)
            str2 = varname * "[" * string(i-1) * "]";
            mystr = mystr + "\t\tint " * str1 * " = " + str2 * ";\n";
        end
    end    
elseif flg==2
    for i = 1:n
        str1 = varname * string(i);
        if contains(ustr, str1)
            str2 = varname * "[" * string(i-1) * " + i*" * string(n) * "]";
            mystr = mystr * "\t\tT " * str1 * " = " * str2 * ";\n";
        end
    end
elseif flg==3
    for i = 1:n        
        str1 = varname * string(i);
        if contains(ustr, str1)
            str2 = varname * "[" * string(i-1) * " + i*" * string(n) * "]";
            mystr = mystr * "\t\tint " * str1 * " = " * str2 * ";\n";            
        end
    end    
end

return mystr;

end
