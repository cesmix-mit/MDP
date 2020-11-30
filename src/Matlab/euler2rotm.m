function R = euler2rotm(eul)
% eul = (alpha, beta, gamma): Euler's angles
% R: rotation matrix

ct = cos(eul);
st = sin(eul);

n = size(eul,1);
R = zeros(n,3,3);

R(:,1,1) = ct(:,1).*ct(:,3).*ct(:,2) - st(:,1).*st(:,3);
R(:,1,2) = -ct(:,1).*ct(:,2).*st(:,3) - st(:,1).*ct(:,3);
R(:,1,3) = ct(:,1).*st(:,2);
R(:,2,1) = st(:,1).*ct(:,3).*ct(:,2) + ct(:,1).*st(:,3);
R(:,2,2) = -st(:,1).*ct(:,2).*st(:,3) + ct(:,1).*ct(:,3);
R(:,2,3) = st(:,1).*st(:,2);
R(:,3,1) = -st(:,2).*ct(:,3);
R(:,3,2) = st(:,2).*st(:,3);
R(:,3,3) = ct(:,2);
R = squeeze(R);
