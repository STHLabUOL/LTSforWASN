function micPos = randomMicPosNearLocSrc(pos_locSrc, dr_range, roomDim, minDistWall)
    % Draw random position near given position of local source
    %
    % In:
    %   posLocSrc (3,1): Reference position of local source
    %   dr_range (2,1): range of desired distance between mic and source
    %   roomDim (3,1): Room dimensions
    %   minDistWall: min. distance to wall
    % Out:
    %   pos (3,1): Mic Position
    %

    micPos = pos_locSrc;
    isOutsideRoom = @(pos) sum(roomDim-minDistWall-pos < 0) > 0 || sum(pos < minDistWall) > 0;

    while isequal(micPos, pos_locSrc) || isOutsideRoom(micPos)
        % cos(theta) and r^3 [instead of theta and r] need to be uniformly distributed for
        % a uniformly random position within the sphere(-ring)
        cos_theta = unifrnd(-1, 1);
        phi = unifrnd(0, 2*pi);
        r3 = unifrnd(dr_range(1).^3, dr_range(2).^3);
        [p1, p2, p3] = sph2cart(phi, acos(cos_theta)-pi/2, r3.^(1/3));    
        micPos = pos_locSrc+[p1, p2, p3];
    end



end