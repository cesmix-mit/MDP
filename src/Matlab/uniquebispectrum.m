function ind = uniquebispectrum(b)
% Find non-zero unique bispectrum components
% b   : 3D array for bispectrum components
% ind : integer array contains indices for non-zero unique bispectrum components

N = size(b,1);
L = N - 1;

% Symmetry conditions for the bispectrum components
% [l2 l1 l]
% [l1 l2 l]
% [l l2 l1]
% [l2 l l1]
% [l1 l l2]
% [l l1 l2]

% use the above symmetry conditions to remove zero and duplicate bispectrum compoments
inc = 0;
for l = 0:L
    for l1 = 0:L
        for l2 = 0:L
            tm = b(l2+1,l1+1,l+1);
            if (abs(tm)>1e-10) && (l2<=l1)
                if (l1==l2)
                    inc = inc + 1;
                elseif (l2<l1) && (l2<l) && (l1<l)
                    inc = inc + 1;                    
                end
            end                
        end
    end
end

ind = zeros(inc,3);
inc = 0;
for l = 0:L
    for l1 = 0:L
        for l2 = 0:L
            tm = b(l2+1,l1+1,l+1);
            if (abs(tm)>1e-10) && (l2<=l1)
                if (l1==l2)
                    inc = inc + 1;
                    ind(inc,:) = [l2 l1 l];
                elseif (l2<l1) && (l2<l) && (l1<l)
                    inc = inc + 1;         
                    ind(inc,:) = [l2 l1 l];
                end
            end                
        end
    end
end

% check 
for l = 0:L
    for l1 = 0:L
        for l2 = 0:L
            tm = b(l2+1,l1+1,l+1);
            im = [l2 l1 l];
            in = ismember(im, ind, 'rows');
            if (abs(tm)>1e-10) && (l2<=l1) && (in==0)
                if (l2~=l1) && (l2~=l) && (l1~=l)
                    [la, lb] = ismember(sort(im), ind, 'rows');
                    if la==1             
                        inb = ind(lb,:);
                        if abs(abs(b(l2+1,l1+1,l+1)/sqrt(2*l+1))-abs(b(inb(1)+1,inb(2)+1,inb(3)+1)/sqrt(2*inb(3)+1))) > 1e-10
                            error("something wrong");                        
                        end
                    else
                        error("something wrong");                        
                    end                    
                elseif (im(1)==im(3))
                    [la, lb] = ismember([im(1) im(3) im(2)], ind, 'rows');
                    if la==1             
                        inb = ind(lb,:);
                        if abs(abs(b(l2+1,l1+1,l+1)/sqrt(2*l+1))-abs(b(inb(1)+1,inb(2)+1,inb(3)+1)/sqrt(2*inb(3)+1))) > 1e-10
                            error("something wrong");                        
                        end
                    else
                        error("something wrong");                        
                    end                    
                elseif (im(2)==im(3))
                    [la, lb] = ismember([im(2) im(3) im(1)], ind, 'rows');
                    if la==1             
                        inb = ind(lb,:);
                        if abs(abs(b(l2+1,l1+1,l+1)/sqrt(2*l+1))-abs(b(inb(1)+1,inb(2)+1,inb(3)+1)/sqrt(2*inb(3)+1))) > 1e-10
                            error("something wrong");                        
                        end
                    else
                        error("something wrong");                        
                    end                    
                else
                    error("something wrong");                        
                end
            end                
        end
    end
end

