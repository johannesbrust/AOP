function [ SM, cidx ] = select_oblique( n, ridx, cidx, col, SM )
%SELECT_OBLIQUE Enumerate all selections of row combinations
%{
    INPUTS:
    n       : Total depth
    ridx    : Current depth in search tree
    cidx    : Column index to write into matrix SM (cidx = 0)
    col     : Selection column to write
    SM      : Selection matrix

    OUTPUTS:
    SM      : Filled Selection matrix

%}

if ridx == n
    % Add columns
    cidx        = cidx +1;
    SM(:,cidx)  = col;

    
%     col(ridx)   = 1;
%     cidx        = cidx +1;
%     SM(:,cidx)  = col;
    
    return;

end

[SM, cidx]      = select_oblique(n,ridx+1, cidx, col,SM);

col_            = col;
col_(ridx+1)    = 1;
[SM, cidx]      = select_oblique(n,ridx+1,cidx, col_,SM);

end


