function [ VARk ] = updateVARk( GAMk, Y, NUk, M)

    VARk = sum( GAMk .* (( Y - NUk ).^2)) / sum( GAMk );
end

