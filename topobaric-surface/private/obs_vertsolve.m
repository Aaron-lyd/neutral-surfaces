function [x, s, t] = obs_vertsolve(SppX, TppX, X, BotK, s, t, x, dfnb, dfnc, s0, t0, tolx) %#codegen
%OBS_VERTSOLVE  Root finding of pressure or depth that matches equation of
%               state with single-valued function.
%
%
% [x, s, t] = obs_vertsolve(SppX, TppX, X, BotK, s, t, x, branch, vafnp, s0, t0, tolx)
% finds the pressure or depth x (within tolerance tolx) and its associated
% practical / Absolute salinity s and potential / Conservative temperature
% t, of a surface on which delta equals that determined by the
% single-valued function determined by dfnb (spline break points) and dfnc
% (spline coefficients), in an ocean whose practical / Absolute salinity
% and potential / Conservative temperature as functions of pressure or
% depth X are given by piecewise polynomials whose coefficients are SppX
% and TppX, and whose knots are X.  The number of valid data points in each
% water column is given by BotK.  The equation of state is given by eos.m
% in the path, taking S, T, and X as its 3 inputs. delta is the in-situ
% density anomaly or specific volume anomaly, defined as eos(s,t,x) -
% eos(s0,t0,x) where s,t are S,T interpolated from X to x, and s0, t0 are
% reference values.  The inputs s and t are not used, but provided so these
% variables may be manipulated in-place.
%
%
% --- Input:
% SppX [O, K-1, N]: coefficients for piecewise polynomial for practical 
%                   / Absolute Salinity in terms of X
% TppX [O, K-1, N]: coefficients for piecewise polynomial for potential 
%                   / Conservative Temperature in terms of X
% X [K, N]: knots for the pressure or depth of the casts
% BotK [1, N]: number of valid data points on each cast
% s [1, N]: initial practical / Absolute salinity on the initial surface
% t [1, N]: initial potential / Conservative temperature on the initial surface
% x [1, N]: initial pressure [dbar] or depth [m] of the surface at each cast
% dfnb [1,B]  : break points for the spline giving delta as a function of x
% dfnc [B-1,D+1]: coefficient matrix for the spline giving delta as a function of x
% s0 [1, 1]: reference S value for delta
% t0 [1, 1]: reference T value for delta
% tolx [1, 1]: tolerance on pressure [dbar] or depth [m] for bisection solver
%
% Note: O is the order of the piecewise polynomials down each cast
%       K is the maximum number of knots in these piecewise polynomials, 
%           i.e. the maximum number of bottles in any cast
%       N is the number of water columns (possibly including land).
%       B is the number of break points in the spline.
%       D is the degree of the spline.
%
% Note: variables can actually be higher dimensional, e.g. N = [ni, nj],
%       and x can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: X must increase along its first dimension.
%
%
% --- Output:
% x [same as input x]: pressure or depth of the updated surface
% s [same as input x]: practical / Absolute salinity of the updated surface
% t [same as input x]: potential / Conservative temperature of the updated surface%
%
%
% --- Acknowledgements:
% The sub-function ppval1 is adapted from MATLAB's function PPVAL.

% --- Copyright:
% This file is part of Neutral Surfaces.
% Copyright (C) 2019  Geoff Stanley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com
% Version   : 2.1.0
%
% Modified by : --
% Date        : --
% Changes     : --

N = numel(x);
Xmat = ~isvector(X);

[K,XN] = size(X);
KX = K * double(XN > 1);
nx = 0;

% Loop over each valid water column
for n = 1:N
    k = BotK(n);
    if k > 1 && isfinite(x(n))
            
        % Select this water column
        SppXn = SppX(:,1:k-1,n);
        TppXn = TppX(:,1:k-1,n);
        if Xmat
          Xn = X(1:k,n);
        else
          Xn = X((1:k).'); % .' is for codegen, so X and (1:k).' both column vectors
        end
        
        lb = Xn(1);
        ub = Xn(k);
        
        % Search for a sign-change, expanding outward from an initial guess 
        [lb, ub] = fzero_guess_to_bounds(@myfcn, x(n), lb, ub, ...
          SppXn, TppXn, Xn, dfnb, dfnc, s0, t0);
        
        if ~isnan(lb)
          % A sign change was discovered, so a root exists in the interval.
          % Solve the nonlinear root-finding problem using Brent's method
          x(n) = fzero_brent(@myfcn, lb, ub, tolx, ...
            SppXn, TppXn, Xn, dfnb, dfnc, s0, t0);
          
          % Interpolate S and T onto the updated surface
          [s(n),t(n)] = ppc_val2(Xn, SppXn, TppXn, x(n));
        else
          x(n) = nan;
          s(n) = nan;
          t(n) = nan;
        end
        
      
    end
    nx = nx + KX;
end

end


function out = myfcn(x, SppX, TppX, X, dfnb, dfnc, s0, t0)
% The difference in delta between the single-valued function and the
% equation of state.
[s,t] = ppc_val2(X, SppX, TppX, x);

out = ppval1(dfnb, dfnc, x) - ( eos(s, t, x) - eos(s0, t0, x) );
end


function v = ppval1(b,c,xx)
% PPVAL1: evaluate 1-dimensional piecewise polynomials.
% b: pp.breaks
% c: pp.coefs
% xx: vector of evaluation locations
% Assumes that d == 1, where d == pp.dim.
% The following code is adapted from MATLAB's PPVAL.

xx = xx(:).';   % Ensure row vector
lx = numel(xx);

% for each evaluation site, compute its breakpoint interval
% (mindful of the possibility that xx might be empty)
if lx > 0
    [~,index] = histc(xx,[-inf,b(2:end-1),inf]); %#ok<HISTC>
else
    index = ones(1,lx);
end

% adjust for NaN. (Inf's are handled naturally by histc).
index(isnan(xx)) = 1;

% now go to local coordinates ...
xx = xx-b(index);

% ... and apply nested multiplication:
v = c(index,1);
for i=2:size(c,2)
    v = xx(:).*v + c(index,i);
end

end
