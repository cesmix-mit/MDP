function [dirname,found] = findexec(filename)

dirname = "";
found = 0;

[status,~] = system("which " + filename);
if status==0    
    dirname = filename;    
    %disp("MDP found " + dirname);
    found = 1;
    return;
end
[status,~] = system("which /usr/bin/" + filename);
if status==0
    dirname = "/usr/bin/" + filename;    
    %disp("MDP found " + dirname);
    found = 1;
    return;
end
[status,~] = system("which /usr/local/bin/" + filename);
if status==0
    dirname = "/usr/local/bin/" + filename;
    %disp("MDP found " + dirname);
    found = 1;
    return;
end
[status,~] = system("which /opt/local/bin/" + filename);
if status==0
    dirname = "/opt/local/bin/" + filename;
    %disp("MDP found " + dirname);
    found = 1;
    return;
end

end

