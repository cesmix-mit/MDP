function mystr = varsassign2(mystr, varname, n, ustr, flg)

if flg==0
    for i = 1:n
        str1 = varname + num2str(i);
        if any(contains(ustr, str1)) 
            str2 = varname + "[" + num2str(i-1) + "]";
            mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
        end
    end
elseif flg==1
    for i = 1:n
        str1 = varname + num2str(i);
        if any(contains(ustr, str1)) 
            str2 = varname + "[" + num2str(i-1) + "]";
            mystr = mystr + "\t\tint " + str1 + " = " + str2 + ";\n";
        end
    end    
elseif flg==2
    for i = 1:n
        str1 = varname + num2str(i);
        if any(contains(ustr, str1)) 
            str2 = varname + "[" + num2str(i-1) + " + i*" + num2str(n) + "]";
            mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";
        end
    end
elseif flg==3
    for i = 1:n        
        str1 = varname + num2str(i);
        if any(contains(ustr, str1)) 
            str2 = varname + "[" + num2str(i-1) + " + i*" + num2str(n) + "]";
            mystr = mystr + "\t\tint " + str1 + " = " + str2 + ";\n";
        end
    end    
end

end

