function OPTS = omega_defaults()
%OMEGA_DEFUALTS  Default options for omega_surface
%
%
% OPTS = omega_defaults()
% returns a struct OPTS containing default options for use in
% omega_surface.

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

OPTS = struct();

%OPTS.AREA_iJ = 1;    % Area [m^2] centred at (I-1/2, J)
%OPTS.AREA_Ij = 1;    % Area [m^2] centred at (I, J-1/2)
OPTS.DIST1_iJ = 1;   % Distance [m] in 1st dimension centred at (I-1/2, J)
OPTS.DIST2_Ij = 1;   % Distance [m] in 2nd dimension centred at (I, J-1/2)
OPTS.DIST2_iJ = 1;   % Distance [m] in 2nd dimension centred at (I-1/2, J)
OPTS.DIST1_Ij = 1;   % Distance [m] in 1st dimension centred at (I, J-1/2)

OPTS.MLX = []; % Do not remove the Mixed Layer

OPTS.FIGS_SHOW = false; % do not show figures

OPTS.INTERPFN = @ppc_linterp; % Use linear interpolation in the vertical dimension.

OPTS.SppX = [];  % Pre-computed interpolation functions.  None given here.
OPTS.TppX = [];  % Pre-computed interpolation functions.  None given here.

OPTS.ITER_MIN = 1;  % minimum number of iterations
OPTS.ITER_MAX = 10; % maximum number of iterations

OPTS.ITER_START_WETTING = 1; % start wetting immediately

OPTS.TOL_LRPD_L1 = 1e-7; % Tolerance in Locally Referenced Potential Density [kg m^-3]

% quit when the L2 change of pressure on surface exceeds this value (set to
% 0 to deactivate), in the same units as X [dbar] or [m].
OPTS.TOL_X_CHANGE_L2 = inf;

% Error tolerance when root-finding to update surface, in the same units as
% X [dbar] or [m].
OPTS.TOL_X_UPDATE = 1e-4;

OPTS.VERBOSE = 1; % show a moderate level of information

OPTS.INTEGRATING_FACTOR = [];  % No integrating factor

OPTS.FILE_ID = 1; % standard output to MATLAB terminal

OPTS.REF_IJ = []; % No reference water column at which the surface is pinned.  Instead, maintain the mean density. 

OPTS.POISSON = true; % Whether to solve the square, symmetric Poisson matrix problem, or a rectangular gradient problem