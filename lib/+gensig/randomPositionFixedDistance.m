function [pos_new] = randomPositionFixedDistance(Room, pos_ref, d_des, d_wall_min)

% Generate a random position for a second mic2  
% with fixed desired distance to reference position of mic1 within a given
% room.
%
% Note that a random position for mic2 in the desired distance for every possible
% mic1 reference-position is guaranteed to exist only if:
% d_des <= 0.5*sqrt(x^2+y^2+z^2), with x,y,z being the coordinates of the
% virtual placement room (given by parRoom.Pos_Reg). 
% If this equation is not fulfilled, there are mic1-positions for
% which no mic2-position in the desired distance exists within the room.
%
% :param Room (1x3): room measurements [x, y, m] in meters
% :param pos_ref (1x3): coordinates of reference position
% :param d_des (scalar): desired distance [m] between positions
% :param d_wall_min (scalar): desired minimum distance from wall

    % Force row vectors
    Room = Room(:).';
    pos_ref = pos_ref(:).';
    
    % Prepare output
    pos_new = zeros(1 ,3);

    % Transform Room and position according to d_wall_min
    Room = Room - d_wall_min*2*ones(1, 3);
    pos_ref = pos_ref - d_wall_min*ones(1, 3);
    
    % Define room-specific max. reachable distance for every dimension
    dx_max_room = Room(1)/2 + abs(Room(1)/2 - pos_ref(1));
    dy_max_room = Room(2)/2 + abs(Room(2)/2 - pos_ref(2));
    dz_max_room = Room(3)/2 + abs(Room(3)/2 - pos_ref(3));

    % Check if desired distance can be reached in the room
    if d_des > sqrt(dx_max_room^2+dy_max_room^2+dz_max_room^2)
        disp(sqrt(dx_max_room^2+dy_max_room^2+dz_max_room^2));
        error('For the given room geometry and reference position, no position in the desired distance exists!');
    end

    % Set X coordinate    
    % -- set min. distance that guaruantees d_des is achievable in 3d room
    dx_min = sqrt(max(0, d_des^2-dy_max_room^2-dz_max_room^2));
    % -- set max. distance either by d_des or room-restricted max. distance
    dx_max = min(d_des, dx_max_room);
    % -- select a random distance from all possible values
    dx = dx_min + rand(1)*(dx_max-dx_min);
    % -- go either positive or negative direction
    if pos_ref(1)+dx > Room(1)
        pos_new(1) = pos_ref(1)-dx;
    elseif pos_ref(1)-dx < 0
        pos_new(1) = pos_ref(1)+dx;
    else
        pos_new(1) = pos_ref(1) + (1-2*round(rand(1)))*dx; %random direction
    end
    
    % Set Y Coordinate
    % -- set min. distance that guaruantees d_des is achievable in 3d room
    dy_min = sqrt(max(0, d_des^2-dx^2-dz_max_room^2));
    % -- set max. distance either by d_des or room-restricted max. distance
    dy_max = min(sqrt(d_des^2-dx^2), dy_max_room);
    % -- select a random distance from all possible values
    dy = dy_min + rand(1)*(dy_max-dy_min);
    % -- go either positive or negative direction
    if pos_ref(2)+dy > Room(2)
        pos_new(2) = pos_ref(2)-dy;
    elseif pos_ref(2)-dy < 0
        pos_new(2) = pos_ref(2)+dy;
    else
        pos_new(2) = pos_ref(2) + (1-2*round(rand(1)))*dy; %random
    end
    
    % Set Z Coordinate
    % -- set exact required z-distance to reach d_des
    dz = sqrt(d_des^2-dx^2-dy^2);
    % -- go either positive or negative direction
    if pos_ref(3)+dz > Room(3)
        pos_new(3) = pos_ref(3)-dz;
    elseif pos_ref(3)-dz < 0
        pos_new(3) = pos_ref(3)+dz;
    else
        pos_new(3) = pos_ref(3) + (1-2*round(rand(1)))*dz; %random
    end
    
    % Transform back to reference room system
    pos_new = pos_new + d_wall_min*ones(1, 3);

end