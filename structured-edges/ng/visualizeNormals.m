function out = visualizeNormals(Z)

if size(Z,3) == 3
  N = Z;
else
  N = getNormals(Z);
end

% hue = (atan2(N(:,:,1), N(:,:,2)) + pi)/(2*pi);
% sat = sqrt(N(:,:,1).^2 + N(:,:,2).^2);
% val = N(:,:,3)/2+0.5;
% V = hsv2rgb(cat(3, hue, sat, val));

N = N(:,:,[3,1,2]);
N(:,:,2:3) = N(:,:,2:3) / 1.25;
V = max(0, min(1, yuv2rgb_simple(N)));

V(isnan(N)) = 0;

if nargout == 0
  imagesc(V);
  imtight;
else
  out = V;
end
