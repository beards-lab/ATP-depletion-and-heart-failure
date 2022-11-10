function arr = insertAt(arr, arr_ins, positions)
    if length(arr_ins) ~= length(positions)
        error("Array to insert must be of the same length");
    end
    arr(positions) = arr_ins;
end