function [ NUk ] = updateNUk( GAMk, Y, M)

    NUk = sum((GAMk .* Y) ) / sum((GAMk));
end

