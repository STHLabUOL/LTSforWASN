function [pos_new] = randomPositionForFixedDistance(room, r0, d_wall_min)

    % This is a sibling function to "randomPositionFixedDistance". For given
    % room measurements, it generates a random position within that room for 
    % which at least one other random position exists in the desired distance d_des.
    %
    % This avoids the problem of needing to draw the first position
    % repeatedly until another (random) position in the desired distance
    % actually exists. (Especially problematic for large desired distances)
    %
    % Input:
    % room (1x3): [x, y, z] room measurements in m
    % r0 (double): desired distance for the future second position
    % d_wall_min (double): min. distance from walls 
    
    % Measurements of virtual placement room
    room = room(:).';
    room = room - 2*d_wall_min*ones(1, 3);

    % Note: Random coordinates are found only for the lower left octant of
    % the room and finally randomly assigned to one of the octants. 
    
    % Draw a
    a_max_0 = room(1)/2;
    a_max = a_max_0;
    is_restricted = r0 >= sqrt(room(2)^2+room(3)^2);
    if is_restricted
        a_max = min(a_max_0, room(1) - sqrt(r0^2-room(2)^2-room(3)^2));
        if a_max < 0
            error('Upper bound for a is negative. No valid results exist for this r0 and room (?)');
        end
    end
    a = unifrnd(0, a_max);
    
    % Draw b    
    b_max_0 = room(2)/2;
    b_max = b_max_0;
    is_restricted = r0 >= sqrt((room(1)-a)^2 + room(3)^2);
    if is_restricted
        b_max = min(b_max_0, room(2) - sqrt(r0^2 - (room(1)-a)^2 - room(3)^2));
        if b_max < 0
            error('Upper bound for b is negative.');
        end
    end
    b = unifrnd(0, b_max);
    
    % Draw c    
    c_max_0 = room(3)/2;
    c_max = c_max_0;
    is_restricted = r0 >= sqrt((room(1)-a)^2 + (room(2)-b)^2);
    if is_restricted
        c_max = min(c_max_0, room(3) - sqrt(r0^2 - (room(1)-a)^2 - (room(2)-b)^2));
    end
    c = unifrnd(0, c_max);
    
    % Select random octant
    flip = rand(1, 3) > 0.5;
    if flip(1)
        a = room(1) - a;
    end
    if flip(2)
        b = room(2) - b;
    end
    if flip(3)
        c = room(3) - c;
    end

    % Transform coordinates from virtual to reference room  
    pos_new = [a, b, c] + d_wall_min*ones(1, 3);

end