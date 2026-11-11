function arr = insertAt(arr, arr_ins, positions)
% Inserts array in given positions. Used for anonymous function to
% subselect set of params. E.g. 
%
% optimfun = @(optimizedParams)costFunc(insertAt(paramSetAll, optimizedParams, selectedIndexes));
% 
    if length(arr_ins) ~= length(positions)
        error("Array to insert must be of the same length");
    end
    arr(positions) = arr_ins;
end