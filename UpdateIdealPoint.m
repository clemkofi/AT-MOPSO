% 
% 
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints, "
% in IEEE Transactions on Evolutionary Computation, 
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Paper URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

function zmin = UpdateIdealPoint(pop, prev_zmin)
    
    if ~exist('prev_zmin', 'var') || isempty(prev_zmin)
        prev_zmin = inf(size(pop(1).Cost));
    end
    
    zmin = prev_zmin;
    for i = 1:numel(pop)
        zmin = min(zmin, pop(i).Cost);
    end

end