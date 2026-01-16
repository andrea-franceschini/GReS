function fractureBasisFromNormal(normal)

%FRACTUREBASISFROMNORMAL
%   Input:
%     normal    : 1x3 vector
%   Output:
%     normalOut : normalized normal
%     lengthVec : unit in-plane vector
%     widthVec  : unit in-plane vector (right-handed)

tol = 1e-14;

% normalize normal
n = normal / norm(normal);

% pick reference axis least aligned with normal
[~,i] = min(abs(n));
a = zeros(1,3);
a(i) = 1;

% first in-plane direction
l = a - dot(a, n) * n;
l = l / norm(l);

% second in-plane direction (right-handed)
w = cross(n, l);
w = w / norm(w);

% safety orientation check
if dot(cross(l, w), n) < 0
    w = -w;
end

fmt = @(v) sprintf('%.15g, %.15g, %.15g', v);
fprintf('normal    = "%s"\n', fmt(n));
fprintf('lengthVec = "%s"\n', fmt(l));
fprintf('widthVec  = "%s"\n', fmt(w));
end

