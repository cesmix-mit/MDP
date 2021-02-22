function [b, c, d, e, m] = uniquesort(a)
% a  (input): an array of n integer elements
% b (putput): unique elements of the original array a
% c (output): the number of counts for each element of the output array b
% d (output): sorted indices of the original array a
% e (output): sorted elements of the original array a
% m (output): the length of the output array b
% Example: a = [9  0  6  2  2  3  1  5  9  3  5  6  7  4  8  9  7  2  7  1]
%          e = [0  1  1  2  2  2  3  3  4  5  5  6  6  7  7  7  8  9  9  9]
%          d = [1  6  19  3  4  17  5  9  13  7  10  2  11  12  16  18  14  0  8  15]+1
%          b = [0  1  2  3  4  5  6  7  8  9]
%          c = [1  2  3  2  1  2  2  3  1  3] -> c = [0  1  3  6  8  9  11  13  16  17  20]        
%          m = 10        
    
% sort array a
[e, d] = sort(a);

[b, c, m] = uniqueelements(e);

function [b, c, m] = uniqueelements(e)  

n = length(e);
t = zeros(1,n);
p = zeros(1,n);
for i = 1:n
    if i<n
        if (e(i+1) - e(i))>0
            t(i) = i;
            p(i) = 1;
        else
            t(i) = 0;
            p(i) = 0;
        end   
    end
    
    if (i==n) 
        t(i) = n;
        p(i) = 1;
    end    
end

p = cumsum(p);    
m = p(n);

b = 0*e;
c = 0*e;
for i = 1:n
    if (t(i) > 0)
        j = p(i);
        k = t(i);
        c(j+1) = k;
        b(j) = e(k);
    end
end


