function unit_vector = DirectionOf(vector)
    if(norm(vector) > 0)
        unit_vector = vector / norm(vector);
    else
        unit_vector = [0; 0];
    end
end

