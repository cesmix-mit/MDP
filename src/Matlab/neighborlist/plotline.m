function plotline(v1, v2)

if length(v1)==1
    plot([v1 0], [v2 0], 'k-', 'LineWidth', 1.5);
elseif length(v1)==2
    plot([v1(1) v2(1)], [v1(2) v2(2)], 'k-', 'LineWidth', 1.5);
elseif length(v1)==3
    plot3([v1(1) v2(1)], [v1(2) v2(2)], [v1(3) v2(3)], 'k-', 'LineWidth', 1.5);    
end