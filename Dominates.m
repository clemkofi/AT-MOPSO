%
% Function to check if a particular particular or member of a population is
% dominated by another 
%

function b = Dominates(x, y)

    if isstruct(x)
        x = x.Cost;
    end
    
    if isstruct(y)
        y = y.Cost;
    end

    b = all(x <= y) && any(x<y);

end