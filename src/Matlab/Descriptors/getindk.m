function indk = getindk(K)

indk = zeros(K*(K+1)/2,2);
count = 0;
for k1 = 1:K
    for k2 = 1:k1      
        count = count + 1;
        indk(count,:) = [k2 k1];
    end
end
