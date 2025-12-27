function LS = LSreinit( X, T, LS )
    % Level Set reinitialization by brute force
    % Allows for multiple level sets.
    % The level sets in the enriched elements are not modified.
    % INPUT
    %   X     nodal coords
    %   T     conectivity matrix
    %   LS    nodal level set 
    %         size(LS) = [number of level sets, number of nodes]
    % OUTPUT
    %   LS    rectified level sets
    
    ngeom = 4;
    for ils = 1:size( LS, 1 )
       fixedNodes = [];
       for ielem = 1:size( T, 1 )
           elem_crossed_LAB = crossedByLevelSet( LS(ils,T(ielem,1:ngeom))' ); 
           if elem_crossed_LAB == true
               fixedNodes = [fixedNodes T(ielem,1:ngeom)];
           end 
       end
       fixedNodes = unique( fixedNodes );
       fixedLS = LS(ils,fixedNodes);
       %
       for inode = 1:size( LS, 2 )
          if isempty( find( fixedNodes == inode, 1 ) ) 
             [dist,ix] = closeNode( X(inode,:), X(fixedNodes,:) );
             s = sign( LS(ils,inode) );
             LS(ils,inode) = s*(sqrt(dist) + abs( fixedLS(ix) ));
          end
       end
    end
end 

function [d,ix] = closeNode( p, v )
    % find the closest point to p in v
    [d,ix] = min( sum( (v - repmat( p, size(v,1), 1 ))'.^2 ) );
end
function [bool,numberOfCrossingLevelSets,firstCrossingLevelSet] = crossedByLevelSet( LS, tol )
    % returns if the element has to be enriched or not
    % INPUT: %   LS   nodal level set values for one element if several columns are present, 
    %                       each one is interpreted as an different level set
    %                   tol  tolerance
    % OUTPUT:    bool (true for elements to be enriched; false: std element)
    %                   numberOfCrossingLevelSets  
   %                    firstCrossingLevelSet      
    if nargin < 2
       tol = 0;
    end
    % All negative level sets
    N = all( LS < tol );
    % All positive level sets
    P = all( LS >= -tol );
    % Split?
    S = ~N & ~P ;
    % The first crossing level set 
    C = find( S == 1, 1 );
    if isempty( C )
       C = 0;
    end
    % M = k, use material k
    % M = 0, use several materials or use the last material
    M = find( P == 1, 1 );
    if isempty( M )
       M = 0;
    end
    bool = (M == 0 & C > 0) | (C < M & C > 0);
    numberOfCrossingLevelSets = length( find( S ) );
    firstCrossingLevelSet = C;
end 
% function booli1 = crossedByLevelSet(LSe)
%     if abs(sum(sign(LSe))) == length(LSe) 
%         booli1 = false; 
%     else
%         booli1 = true; 
%     end 
% end 