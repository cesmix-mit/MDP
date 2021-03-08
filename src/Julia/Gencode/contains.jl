function contains(ustr, s)

itsin = false;    
if occursin(s, string(ustr))
    itsin = true;
end

return itsin;

end

