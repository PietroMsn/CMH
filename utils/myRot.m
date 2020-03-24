function out = myRot(dato3d, angles)

    if not(angles(1)==0)
        dato3d(:,2) = dato3d(:,2)*cos(angles(1)) - dato3d(:,3)*sin(angles(1));
        dato3d(:,3) = dato3d(:,2)*sin(angles(1)) + dato3d(:,3)*cos(angles(1));
    end
    if not(angles(2)==0)
        dato3d(:,1) = dato3d(:,1)*cos(angles(2)) + dato3d(:,3)*sin(angles(2));
        dato3d(:,3) = dato3d(:,3)*cos(angles(2)) - dato3d(:,1)*sin(angles(2));
    end
    if not(angles(3)==0)
        dato3d(:,1) = dato3d(:,1)*cos(angles(3)) - dato3d(:,2)*sin(angles(3));
        dato3d(:,2) = dato3d(:,1)*sin(angles(3)) + dato3d(:,2)*cos(angles(3));
    end

    out = dato3d;
end