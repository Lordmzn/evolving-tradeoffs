function u = Test_SOP(pol, level)

A = pol.theta(1);
B = pol.theta(2);

if length(level) == 1
    if level <=0
        u = 0;
        return
    elseif level <= A
        u = level;
        return
    elseif level <= B
        u = max(A, level - 100);
        return
    else
        u = level;
    end
else
    
    u(level <= 0) = 0;
    
    u(level > 0 & level <= A) = level(level > 0 & level <= A);
    
    u(level > A & level <= B) = max(repmat(A, size(level(level > A & level <= B))), ...
        level(level > A & level <= B) - 100);
    
    u(level > B) = level(level > B);
    
end