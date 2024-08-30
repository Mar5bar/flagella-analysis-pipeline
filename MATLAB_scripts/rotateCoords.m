function [rx,ry] = rotateCoords(x,y,theta)
% Rotate coords in row vectors x,y by theta.
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rotated = R*[x;y];
    rx = rotated(1,:);
    ry = rotated(2,:);
end