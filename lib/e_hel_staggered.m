function [e_hel, e_hel_rms, e_hel_u, e_hel_v, e_hel_u_iph, e_hel_v_jph] = e_hel_staggered(x, X, u, v, dx, dy, SppX, TppX)


[ni,nj] = size(x);
nk = size(X,1);

% Just-in-Time code generation
ntp_midpoint_to_casts_codegen(nk, ni, nj, 1);  % Currently only works for X a vector

tolx = 1e-5;
% Solve nonlinear problem for each pair of adjacent casts
[dxi, dxj] = ntp_midpoint_to_casts_mex(SppX, TppX, X, tolx, x);
% [dxi, dxj] = ntp_midpoint_to_casts(SppX, TppX, X, tolx, x);


im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west

ip1 = @(D) circshift(D, [-1 0]);
jp1 = @(D) circshift(D, [0 -1]);

lead1 = @(x) reshape(x, [1 size(x)]);

% slope error = grad_n z - grad_a z,  where z < 0
% Negative sign added here because we've been x = -z > 0 people. 
sx = -(dxi - (x - im1(x)) ) ./ dx;
sy = -(dxj - (x - jm1(x)) ) ./ dy;

% z, p, u, v on the u and v grid

x_on_U = (x + im1(x)) / 2;
x_on_V = (x + jm1(x)) / 2;

u_on_U = interp1qn(lead1(x_on_U), X, u);
v_on_V = interp1qn(lead1(x_on_V), X, v);

% Calculate e^hel on (i-1/2, j):
e_hel_u = u_on_U.*sx;
% Calculate e^hel on (i, j-1/2):
e_hel_v = v_on_V.*sy;
% Calculate e^hel on (i+1/2, j):
e_hel_u_iph = ip1(e_hel_u);
% Calculate e^hel on (i, j+1/2):
e_hel_v_jph = jp1(e_hel_v);

e_hel = (e_hel_u + e_hel_v + e_hel_u_iph + e_hel_v_jph)/2;

e_hel_rms = e_hel - nanmean(e_hel(:));
good = ~isnan(e_hel_rms);
e_hel_rms = rms(e_hel_rms(good));

end



