function out = iga_thb_vms_john53_v7 (varargin)
%IGA_THB_VMS_JOHN53_V7  Adaptive THB-spline IGA with VMS estimator for
%   V. John, CMAME 190 (2000) Ex. 5.3 (circular inner layer).
%     -eps*Lap(u) + b.grad(u) + c*u = f in (0,1)^2,  u=0 on dOmega
%     b=(2,3)^T, c=2, eps=1e-4;  u = 16 x(1-x) y(1-y)*(1/2+atan(A)/pi)
%   Runs UNIFORM then ADAPTIVE (THB, Doerfler) per degree in prm.degrees.
%   eta_K = tau_K ||R||_{L2(K)},  tau_K = min(ell/(sqrt(2)|b|), ell^2/(3 sqrt(10) eps))
%   Two drivers: 'vms' marks on eta_K; 'gradind' marks on ||grad u_h||_K.
%   USAGE
%   Requires GeoPDEs 3.x (hierarchical/adaptivity) and the NURBS toolbox.

%%                          USER PARAMETERS
prm.eps        = 1e-4;
prm.b          = [2; 3];
prm.c          = 2;
prm.r0         = 0.25;
prm.xc         = [0.5; 0.5];
prm.degrees    = [1 2 3 4 5];    % -> p = 1..5,  C^{p-1}
prm.nsub0      = 4;
prm.niter      = [];
prm.niter_uniform  = 6;
prm.niter_adaptive = 40;
prm.max_ndof   = 2000;
prm.refinement = 'adaptive';     % 'uniform' | 'adaptive'
prm.test_op    = 'sgs';          % 'sgs' (= MS, the paper) | 'supg' | 'gls'
prm.stabilise  = true;
prm.drivers    = {'vms','gradind'};
% --- hierarchical space ----------------------------------------------
prm.space_type = 'standard';     % 'standard' 
prm.truncated  = true;           % THB, 
% --- estimator, ----------------------------
prm.est_primary = 'local';       %  is THE estimator; I_eff = eta/||e||
prm.norm_type   = 'hs';          % 'hs' = paper (4.7),(4.8) | 'operator'
prm.gamma_safe  = 1.00;          % KEEP AT 1: no manufactured reliability
prm.compute_tr  = false;         % characteristic (4.5) variant, off
prm.bl_correction = false;
prm.est_cert    = true;
prm.chord       = 'stream_harm'; % reduction length ell of (4.10):
                                 %   'max'         = eq (4.10) as published
                                 %   'stream_harm' = set here
                                 %   'flux_rms'    = eq (4.12)
prm.tau_dif     = 'chord';       % 'chord' | 'geometric' | 'min_edge'
% --- diagnostics  --------------
prm.quad_check = false;
prm.diag_proj  = false;
% --- streamline transport solver -----------
prm.tr_Ktr     = 4;
prm.tr_NLmin   = 512;
prm.tr_NLmax   = 8000;
prm.tr_nq      = 6;
prm.tr_dsfac   = 0.5;
prm.tr_dscoarse= 0.05;
prm.tr_band    = 6;
prm.tr_block   = 256;
prm.tr_check   = false;
% --- marking ----------------------------------------------------------
prm.mark_strategy = 'dorfler';
prm.mark_object= 'functions';
prm.mark_on    = 'local';        % the paper's (4.15) marks on eta_e of (4.13)
prm.theta      = 0.40;
prm.adm_class  = 2;
prm.adm_q      = [];
% --- stabilisation ----------------------------------------------------
prm.h_stab     = 'max_chord';
% --- gradient indicator (4.16) ---------------------------------------
prm.gradind_scale = 'none';
% --- quadrature -------------------------------------------------------
prm.nquad      = [];             % [] -> p+2 or can be change to 4,5
prm.quad_mode  = 'subcell';
prm.lay_bandw  = 4;
prm.lay_nsub   = 16;
prm.lay_nq     = 6;
prm.loc_tol    = 1e-3;
% --- linear solver ----------------------------------------------------
prm.iter_refine = 2;
prm.do_condest  = true;
prm.condest_max = 15000;
prm.err_floor_factor = 3;
prm.stop_at_floor    = false;
% --- memory -----------------------------------------------------------
prm.mem_budget = 120e6;          % bytes per element-evaluation chunk
% --- experiments ------------------------------------------------------
prm.damkohler  = [];
prm.selftest   = false;
prm.verbose    = true;
% --- figures ----------------------------------------------------------
prm.plot_npts   = 260;
prm.plots.error    = true;       % L2 error vs dof, all runs, one axes
prm.plots.ieff     = true;       % I_eff   vs dof, all runs, one axes
prm.plots.maps     = true;       % local effectivity map, ONE FIGURE PER
                                 % (degree, driver) -- see plot_local_maps
prm.plots.mesh     = true;       % final adaptive meshes, one grid figure
prm.plots.solution = true;       % surface plots, ONE FIGURE PER
                                 % (degree, driver) -- see plot_solution_surfs
prm.fig_export  = '';            % folder for pdf/png output; '' = do not save
prm.fig_fmt     = {'pdf','png'};
prm.tex_export  = '';            % file for booktabs versions of the tables,
                                 % e.g. 'tables.tex'; '' = do not write
prm.fig_font    = 'Times New Roman';
prm.fig_fsize   = 11;
prm.xpad        = 2.2;           % padding factor of the logarithmic dof axis

prm = parse_options (prm, varargin);
prm = validate_options (prm);

if prm.selftest
  out = run_selftest (prm);
  if nargout == 0, clear out; end
  return;
end
if ~isempty (prm.damkohler)
  out = run_damkohler (prm);
  if nargout == 0, clear out; end
  return;
end

drivers = active_drivers (prm);
R = [];  pb = [];  niter = [];
for j = 1:numel (drivers)
  [Rj, pb, niter] = run_sweep (prm, drivers{j});
  if isempty (R), R = Rj; else, R = [R, Rj]; end   %#ok<AGROW>
end
SUM = build_summary (R, prm);
out = struct ('runs', {R}, 'summary', {SUM}, 'prm', prm, 'niter', niter, ...
              'drivers', {drivers});
print_tables (R, prm);
if prm.plots.error,    plot_error_all    (R, prm); end
if prm.plots.ieff,     plot_ieff_all     (R, prm); end
if prm.plots.maps,     plot_local_maps   (R, prm); end
if prm.plots.mesh,     plot_final_meshes (R, prm); end
if prm.plots.solution, plot_solution_surfs (R, prm, pb); end
if nargout == 0, clear out; end
end

%%  CHUNKED LEVEL EVALUATION

function nc = chunk_size (hspace, ilev, nqn, prm)
%CHUNK_SIZE  Elements per chunk so hessian storage stays within prm.mem_budget.
sp  = hspace.space_of_level(ilev);
nsh = prod (double (sp.degree) + 1);
bpe = 8 * nqn * nsh * 8 * 3;
nc  = max (1, floor (prm.mem_budget / max (bpe, 1)));
end

function [msh_c, sp_c] = eval_elems (hmsh, hspace, ilev, elems)
%EVAL_ELEMS  Evaluate shape functions/gradients/laplacians on a subset of active elements at ilev.
msh_c = msh_evaluate_element_list (mesh_of_level (hmsh, ilev), elems);
sp_c  = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_c, ...
                  'value', true, 'gradient', true, 'laplacian', true);
if ~isfield (sp_c, 'shape_function_laplacians')
  error ('iga:lap', ['GeoPDEs did not return shape_function_laplacians at ' ...
    'p = %d.  The residual would silently lose its -eps*Lap term.'], ...
    hspace.space_of_level(ilev).degree(1));
end
sp_c = remap_to_csub (sp_c, hmsh, hspace, ilev);
end

function sp_c = eval_conn (hmsh, hspace, ilev, elems)
%EVAL_CONN  Connectivity only, for MARK_FUNCTIONS.
msh_c = msh_evaluate_element_list (mesh_of_level (hmsh, ilev), elems);
sp_c  = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_c, ...
                  'value', false, 'gradient', false);
sp_c = remap_to_csub (sp_c, hmsh, hspace, ilev);
end

function sp_c = remap_to_csub (sp_c, hmsh, hspace, ilev)
nrows = size (hspace.Csub{ilev}, 1);
if nrows == sp_c.ndof, return; end
done = false;
if exist ('change_connectivity_localized_Csub', 'file')
  try
    sp_c = change_connectivity_localized_Csub (sp_c, hspace, ilev);
    done = (sp_c.ndof == nrows) && all (sp_c.connectivity(:) > 0);
  catch
    done = false;
  end
end
if ~done
  ridx = csub_row_indices (hspace, ilev, nrows);
  if isempty (ridx), ridx = row_map_from_connectivity (hmsh, hspace, ilev, nrows); end
  if isempty (ridx)
    error ('iga:csub','Could not map level-%d connectivity into Csub{%d}.', ilev, ilev);
  end
  [tf, loc] = ismember (sp_c.connectivity, ridx);
  if ~all (tf(:))
    error ('iga:csub','Could not map level-%d connectivity into Csub{%d}.', ilev, ilev);
  end
  sp_c.connectivity = loc;
  sp_c.ndof = nrows;
end
end


%%  LOCAL TENSOR EVALUATION ON ONE KNOT SPAN

function [D, cols] = span_ders (knots, p, x, nd)
%SPAN_DERS  Derivatives 0..nd of the p+1 active B-splines at points X (all in one span).
knots = knots(:).';
n = numel (knots) - p - 1;
x = min (max (x(:).', knots(1)), knots(end));
nx = numel (x);
s  = findspan (n-1, p, x, knots);
nu = min (nd, p);
Dd = basisfunder (s, p, x, knots, nu);
D  = cell (1, nd+1);
for d = 0:nd
  if d > nu, D{d+1} = zeros (nx, p+1);
  else,      D{d+1} = reshape (Dd(:, d+1, :), nx, p+1);
  end
end
c0 = double (s(1)) - p + 1;
if any (double (s) ~= double (s(1)))
  error ('iga:span','span_ders was given points from more than one knot span.');
end
cols = min (max (c0:(c0+p), 1), n);
end


%% ESTIMATOR SUPPORT: hierarchical location and scattered evaluation
function d = empty_dg ()
d = struct ('NL',0,'nseg',0,'area',0,'ppe',0,'hmin',0,'nstep',0,'ndrop',0);
end

function M = cumgauss (nq)
%CUMGAUSS  M(i,j) = int_{-1}^{xi_i} L_j(t) dt; cached by nq.
persistent CACHE
if isempty (CACHE), CACHE = containers.Map ('KeyType','double','ValueType','any'); end
if CACHE.isKey (nq), M = CACHE(nq); return; end
[xg, ~] = gauss_ref (nq);
V = zeros (nq);
for m = 1:nq, V(:,m) = xg(:).^(m-1); end
Li = inv (V);
M = zeros (nq);
for i = 1:nq
  for j = 1:nq
    t = 0;
    for m = 1:nq, t = t + Li(m,j)*(xg(i)^m - (-1)^m)/m; end
    M(i,j) = t;
  end
end
if max (abs (sum (M,2).' - (xg + 1))) > 1e-10
  error ('iga:cumgauss','cumulative Gauss matrix failed its identity test.');
end
CACHE(nq) = M;
end

function LUT = build_level_lut (hmsh, hspace)
nl = hmsh.nlevels;
LUT = repmat (struct ('bx',[],'by',[],'nd',[0 0],'pos',[],'off',0, ...
                      'nel',0,'hmin',inf,'uni',false,'x0',0,'hx',0, ...
                      'y0',0,'hy',0,'knx',[],'kny',[],'dgx',0,'dgy',0, ...
                      'Cm',[]), 1, nl);
off = 0;
for l = 1:nl
  brk = level_breaks (hmsh, hspace, l);
  bx = unique (brk{1}(:));  by = unique (brk{2}(:));
  nd = [numel(bx)-1, numel(by)-1];
  pos = zeros (prod (nd), 1);
  nel = numel (hmsh.active{l});
  if nel > 0, pos(hmsh.active{l}(:)) = (1:nel).'; end
  LUT(l).bx = bx;  LUT(l).by = by;  LUT(l).nd = nd;
  LUT(l).pos = pos;  LUT(l).off = off;  LUT(l).nel = nel;
  dx = diff (bx);  dy = diff (by);
  LUT(l).uni = (max (abs (dx - dx(1))) < 1e-13*dx(1)) && ...
               (max (abs (dy - dy(1))) < 1e-13*dy(1));
  LUT(l).x0 = bx(1);  LUT(l).hx = dx(1);
  LUT(l).y0 = by(1);  LUT(l).hy = dy(1);
  if nel > 0, LUT(l).hmin = min (min (dx), min (dy)); else, LUT(l).hmin = inf; end
  sp = hspace.space_of_level(l);
  LUT(l).knx = sp.knots{1};  LUT(l).kny = sp.knots{2};
  LUT(l).dgx = sp.degree(1); LUT(l).dgy = sp.degree(2);
  off = off + nel;
end
end

function LUT = attach_coefficients (LUT, hmsh, hspace, u)
ndofs = 0;
for l = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(l);
  if LUT(l).nel == 0, LUT(l).Cm = []; continue; end
  sp = hspace.space_of_level(l);
  cl = hspace.Csub{l} * u(1:ndofs);
  cf = level_coefficients (hmsh, hspace, l, cl, sp.ndof);
  n1 = numel (LUT(l).knx) - LUT(l).dgx - 1;
  n2 = numel (LUT(l).kny) - LUT(l).dgy - 1;
  LUT(l).Cm = reshape (cf, n1, n2);
end
end

function [gl, xlo, xhi, ylo, yhi] = locate_hier_batch (LUT, X, Y)
%LOCATE_HIER_BATCH  Global (level-major) active-element index and bounding box for each point.
M = numel (X);
gl = zeros (M,1);  xlo = gl; xhi = gl; ylo = gl; yhi = gl;
todo = true (M,1);
for l = numel (LUT):-1:1
  if ~any (todo) || LUT(l).nel == 0, continue; end
  id = find (todo);
  ix = cell_index (LUT(l), X(id), 1);
  iy = cell_index (LUT(l), Y(id), 2);
  lin = ix + (iy-1)*LUT(l).nd(1);
  pp = LUT(l).pos(lin);
  hit = pp > 0;
  if any (hit)
    j = id(hit);
    gl(j)  = LUT(l).off + pp(hit);
    xlo(j) = LUT(l).bx(ix(hit));    xhi(j) = LUT(l).bx(ix(hit)+1);
    ylo(j) = LUT(l).by(iy(hit));    yhi(j) = LUT(l).by(iy(hit)+1);
    todo(j) = false;
  end
end
end

function i = cell_index (L, v, dir)
if dir == 1, b = L.bx; n = L.nd(1); v0 = L.x0; h = L.hx;
else,        b = L.by; n = L.nd(2); v0 = L.y0; h = L.hy; end
if L.uni
  i = floor ((v - v0)/h) + 1;
else
  i = discretize_ (v, b, n);
end
i = min (max (i, 1), n);
i = i(:);
end
function i = discretize_ (v, b, n)
i = sum (bsxfun (@ge, v(:), b(1:end-1).'), 2);
i(v(:) >= b(end)) = n;
end

function R = hier_residual_scattered (LUT, pb, X, Y)
%HIER_RESIDUAL_SCATTERED  R = -eps*Lap(u_h)+b.grad(u_h)+c*u_h-f at arbitrary points.
M = numel (X);
R = zeros (M,1);
gl = locate_hier_batch (LUT, X, Y);
lev = zeros (M,1);
for l = 1:numel (LUT)
  if LUT(l).nel == 0, continue; end
  lev (gl > LUT(l).off & gl <= LUT(l).off + LUT(l).nel) = l;
end
for l = 1:numel (LUT)
  id = find (lev == l);
  if isempty (id) || isempty (LUT(l).Cm), continue; end
  blk = 40000;
  for a = 1:blk:numel (id)
    z = id(a:min (a+blk-1, numel (id)));
    [V0,Vx,Vy,Vxx,Vyy] = tp_eval_scattered (LUT(l), X(z), Y(z));
    R(z) = -pb.eps*(Vxx+Vyy) + pb.b(1)*Vx + pb.b(2)*Vy + pb.c*V0 ...
           - pb.f (X(z), Y(z));
  end
end
end

function [V0,Vx,Vy,Vxx,Vyy] = tp_eval_scattered (L, X, Y)
knx = L.knx;  kny = L.kny;  px = L.dgx;  py = L.dgy;  Cm = L.Cm;
nx = numel (knx) - px - 1;  ny = numel (kny) - py - 1;
X = min (max (X(:), knx(1)), knx(end));
Y = min (max (Y(:), kny(1)), kny(end));
M = numel (X);
nux = min (2, px);  nuy = min (2, py);
sx = findspan (nx-1, px, X.', knx);   Dx = basisfunder (sx, px, X.', knx, nux);
sy = findspan (ny-1, py, Y.', kny);   Dy = basisfunder (sy, py, Y.', kny, nuy);
Bx = cell (1,3);  By = cell (1,3);
for k = 0:2
  if k <= nux, Bx{k+1} = reshape (Dx(:,k+1,:), M, px+1); else, Bx{k+1} = zeros (M,px+1); end
  if k <= nuy, By{k+1} = reshape (Dy(:,k+1,:), M, py+1); else, By{k+1} = zeros (M,py+1); end
end
ix = min (max (bsxfun (@plus, double(sx(:)) - px, 0:px) + 1, 1), nx);
iy = min (max (bsxfun (@plus, double(sy(:)) - py, 0:py) + 1, 1), ny);
V0 = zeros (M,1); Vx = V0; Vy = V0; Vxx = V0; Vyy = V0;
for b = 1:py+1
  cb = Cm (bsxfun (@plus, ix, (iy(:,b)-1)*nx));   % M x (px+1)
  t0 = sum (cb .* Bx{1}, 2);
  t1 = sum (cb .* Bx{2}, 2);
  t2 = sum (cb .* Bx{3}, 2);
  V0  = V0  + t0 .* By{1}(:,b);
  Vx  = Vx  + t1 .* By{1}(:,b);
  Vy  = Vy  + t0 .* By{2}(:,b);
  Vxx = Vxx + t2 .* By{1}(:,b);
  Vyy = Vyy + t0 .* By{3}(:,b);
end
end


%% ======================================================================
%  TAU  (paper (4.7)-(4.10))
%  ======================================================================
function [CA, CD] = vms_constants (prm)
%VMS_CONSTANTS  CA, CD for 'hs' (Props 4.1-4.2) or 'operator' (2/pi, 1/pi^2).
switch lower (prm.norm_type)
  case 'operator', CA = 2/pi;        CD = 1/pi^2;
  case 'hs',       CA = 1/sqrt (2);  CD = 1/(3*sqrt (10));
  otherwise, error ('iga:norm','norm_type must be ''operator'' or ''hs''.');
end
end

function ell = stream_chord (hx, hy, ahat, mode)
a1 = abs (ahat(1));  a2 = abs (ahat(2));  tiny = 1e-14;
Lx = hx / max (a1, tiny);
Ly = hy / max (a2, tiny);
A  = min (Lx, Ly);   B = max (Lx, Ly);
switch lower (mode)
  case 'max',      ell = A;                                  % paper (4.10)
  case 'flux_rms', ell = A * sqrt (max (1 - A/(2*B), 0));    % paper (4.12)
  case 'stream_harm', ell = 2 / (a1/max(hx,tiny) + a2/max(hy,tiny));
  case 'flow_extent', ell = a1*hx + a2*hy;
  case 'mean_chord',  ell = hx*hy / max (a2*hx + a1*hy, tiny);
  case 'geometric',   ell = sqrt (hx*hy);
  otherwise, error ('iga:chord','Unknown chord mode ''%s''.', mode);
end
end

function t = tau_up_elem (hx, hy, pb, prm)
[CA, CD] = vms_constants (prm);
ell  = stream_chord (hx, hy, pb.ahat, prm.chord);
ellm = stream_chord (hx, hy, pb.ahat, 'max');
adv  = max (CA * ell / pb.bnorm, realmin);
switch lower (prm.tau_dif)
  case {'chord','max_chord'}, dif = CD * ellm^2 / pb.eps;
  case 'geometric',           dif = CD * (hx*hy) / pb.eps;
  case 'min_edge',            dif = CD * min (hx,hy)^2 / pb.eps;
  otherwise, error ('iga:taudif','Unknown tau_dif ''%s''.', prm.tau_dif);
end
t = min (adv, max (dif, realmin));
end

function h = elem_length (hx, hy, ahat, mode)
a1 = abs (ahat(1));  a2 = abs (ahat(2));  tiny = 1e-14;
switch lower (mode)
  case 'max_chord',   h = min (hx / max (a1,tiny), hy / max (a2,tiny));
  case 'mean_chord',  h = hx * hy / max (a2*hx + a1*hy, tiny);
  case 'geometric',   h = sqrt (hx * hy);
  case 'max_edge',    h = max (hx, hy);
  case 'min_edge',    h = min (hx, hy);
  otherwise, error ('iga:hmode','Unknown element-length mode ''%s''.', mode);
end
end

function v = xi_coth (a)
if a < 1e-8, v = a/3; else, v = coth (a) - 1/a; end
end

%% ======================================================================
%  ESTIMATE + TRUE ERROR
%  ======================================================================
function E = estimate_and_error (hmsh, hspace, u, pb, prm)
ntot = hmsh.nel;
r2 = zeros (ntot,1);  e2 = zeros (ntot,1);  g2 = zeros (ntot,1);
q2 = zeros (ntot,1);  r2p = zeros (ntot,1); e2p = zeros (ntot,1);
lev_of = zeros (ntot,1);  pos_of = zeros (ntot,1);
cxK = zeros (ntot,1);  cyK = zeros (ntot,1);
hxK = zeros (ntot,1);  hyK = zeros (ntot,1);
sub  = strcmpi (prm.quad_mode, 'subcell');
band = prm.lay_bandw * pb.delta;
off = 0;  ndofs = 0;
for ilev = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(ilev);
  nel_lev = hmsh.nel_per_level(ilev);
  if nel_lev == 0, continue; end
  u_lev = hspace.Csub{ilev} * u(1:ndofs);
  sp_tp = hspace.space_of_level(ilev);
  if sub
    coef  = level_coefficients (hmsh, hspace, ilev, u_lev, sp_tp.ndof);
    knx = sp_tp.knots{1};  kny = sp_tp.knots{2};
    dgx = sp_tp.degree(1); dgy = sp_tp.degree(2);
    n1 = numel (knx) - dgx - 1;  n2 = numel (kny) - dgy - 1;
    Cm = reshape (coef, n1, n2);
  end
  act = hmsh.active{ilev}(:);
  nc  = chunk_size (hspace, ilev, prm.nquad^2, prm);
  for a0 = 1:nc:nel_lev
    sel  = a0 : min (a0+nc-1, nel_lev);
    [msh_c, sp_c] = eval_elems (hmsh, hspace, ilev, act(sel));
    [hxe, hye, cxe, cye] = elem_geom_of_level (hmsh, hspace, ilev, act(sel));
    nsh = sp_c.nsh_max;  nqn = msh_c.nqn;
    for iel = 1:numel (sel)
      k = off + sel(iel);
      w = msh_c.quad_weights(:,iel) .* msh_c.jacdet(:,iel);
      x = reshape (msh_c.geo_map(1,:,iel), nqn, 1);
      y = reshape (msh_c.geo_map(2,:,iel), nqn, 1);
      [N, Gx, Gy, Lp] = shape_blocks (sp_c, iel, nqn, nsh);
      ue = u_lev(sp_c.connectivity(:,iel));
      uh = N*ue;  gx = Gx*ue;  gy = Gy*ue;  lp = Lp*ue;
      Rh = -pb.eps*lp + pb.b(1)*gx + pb.b(2)*gy + pb.c*uh - pb.f (x,y);
      r2p(k) = max (sum (w .* Rh.^2), 0);
      e2p(k) = max (sum (w .* (pb.u (x,y) - uh).^2), 0);
      q2(k)  = max (sum (w .* uh.^2), 0);
      r2(k)  = r2p(k);  e2(k) = e2p(k);
      g2(k)  = max (sum (w .* (gx.^2 + gy.^2)), 0);
      lev_of(k) = ilev;  pos_of(k) = act(sel(iel));
      cxK(k) = cxe(iel);  cyK(k) = cye(iel);
      hxK(k) = hxe(iel);  hyK(k) = hye(iel);
      if sub && circle_crosses (cxe(iel), cye(iel), hxe(iel), hye(iel), pb, band)
        ns = subcell_count (max (hxe(iel), hye(iel)), pb, prm);
        [xq, wxq] = subcell_rule (cxe(iel)-hxe(iel)/2, cxe(iel)+hxe(iel)/2, ns, prm.lay_nq);
        [yq, wyq] = subcell_rule (cye(iel)-hye(iel)/2, cye(iel)+hye(iel)/2, ns, prm.lay_nq);
        [Bx, cix] = span_ders (knx, dgx, xq, 2);
        [By, ciy] = span_ders (kny, dgy, yq, 2);
        Cl = Cm(cix, ciy);                       % (dgx+1) x (dgy+1), LOCAL
        V0  = Bx{1}*Cl*By{1}.';   Vx  = Bx{2}*Cl*By{1}.';
        Vy  = Bx{1}*Cl*By{2}.';   Vxx = Bx{3}*Cl*By{1}.';
        Vyy = Bx{1}*Cl*By{3}.';
        [XG, YG] = ndgrid (xq, yq);
        WG = wxq(:) * wyq(:).';
        RG = -pb.eps*(Vxx+Vyy) + pb.b(1)*Vx + pb.b(2)*Vy + pb.c*V0 - pb.f (XG,YG);
        r2(k) = max (sum (sum (WG .* RG.^2)), 0);
        e2(k) = max (sum (sum (WG .* (pb.u (XG,YG) - V0).^2)), 0);
        q2(k) = max (sum (sum (WG .* V0.^2)), 0);
        g2(k) = max (sum (sum (WG .* (Vx.^2 + Vy.^2))), 0);
      end
    end
    clear msh_c sp_c;
  end
  off = off + nel_lev;
end

tauK = zeros (ntot,1);
for k = 1:ntot, tauK(k) = tau_up_elem (hxK(k), hyK(k), pb, prm); end

E.eta_locK  = tauK .* sqrt (r2);          % the paper's (4.13)
E.err_K     = sqrt (e2);
E.err_Kp    = sqrt (e2p);
E.eta_gradK = gradind_scale_factor (hxK, hyK, prm) .* sqrt (g2);   % (4.16)

if ~prm.compute_tr
  E.eta_trK = zeros (ntot,1);  E.dg = empty_dg ();  E.dg.area = 1;
  E.tr_rich = NaN;  E.dg.NL2 = NaN;
else
  LUT = build_level_lut (hmsh, hspace);
  LUT = attach_coefficients (LUT, hmsh, hspace, u);
  [E.eta_trK, E.dg] = transport_stream (LUT, pb, prm, ntot, []);
  if prm.tr_check
    q2p = prm;  q2p.tr_NLmin = 2*prm.tr_NLmin;  q2p.tr_Ktr = 2*prm.tr_Ktr;
    q2p.tr_nq = prm.tr_nq + 3;  q2p.tr_dsfac = prm.tr_dsfac/2;
    q2p.tr_dscoarse = prm.tr_dscoarse/2;
    [t2, d2] = transport_stream (LUT, pb, q2p, ntot, []);
    E.tr_rich = sqrt (sum (t2.^2));
    E.dg.NL2  = d2.NL;
  else
    E.tr_rich = NaN;  E.dg.NL2 = NaN;
  end
end

switch lower (prm.est_primary)
  case 'transport', E.eta_K = prm.gamma_safe * E.eta_trK;
  case 'local',     E.eta_K = prm.gamma_safe * E.eta_locK;
  case 'max',       E.eta_K = prm.gamma_safe * max (E.eta_trK, E.eta_locK);
  otherwise, error ('iga:prim','Unknown est_primary ''%s''.', prm.est_primary);
end

E.eta      = sqrt (sum (E.eta_K.^2));
E.eta_loc  = sqrt (sum (E.eta_locK.^2));
E.eta_tr   = sqrt (sum (E.eta_trK.^2));
E.eta_grad = sqrt (sum (E.eta_gradK.^2));
E.err      = sqrt (sum (E.err_K.^2));
E.err_p    = sqrt (sum (E.err_Kp.^2));
E.rnorm    = sqrt (sum (r2));
E.unorm    = sqrt (sum (q2));
E.eta_cert = E.rnorm / pb.c0;
E.lev_of = lev_of;  E.pos_of = pos_of;
E.cx = cxK;  E.cy = cyK;  E.hx = hxK;  E.hy = hyK;
E.tau = tauK;
E.rK  = sqrt (r2);        % kept so any tau variant can be scored offline

% ---- where the error and the estimator live: layer vs elsewhere -------
rc = hypot (cxK - pb.xc(1), cyK - pb.xc(2));
isl = abs (rc - pb.r0) < 0.5*hypot (hxK, hyK) + 2*pb.delta;
E.frac_err_lay = sum (e2(isl)) / max (sum (e2), realmin);
E.frac_eta_lay = sum (E.eta_locK(isl).^2) / max (sum (E.eta_locK.^2), realmin);

if prm.quad_check
  [E.eta_ref, E.err_ref] = ref_quad_norms (hmsh, hspace, u, pb, prm, tauK);
else
  E.eta_ref = NaN;  E.err_ref = NaN;
end
end

function [eta_ref, err_ref] = ref_quad_norms (hmsh, hspace, u, pb, prm, tauK)
%REF_QUAD_NORMS  eta and ||e|| with enriched sub-cell quadrature; checks quadrature saturation.
q = prm;
q.lay_nq   = prm.lay_nq + 4;
q.lay_nsub = 2*prm.lay_nsub;
r2 = zeros (hmsh.nel,1);  e2 = zeros (hmsh.nel,1);
off = 0;  ndofs = 0;
for ilev = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(ilev);
  nel_lev = hmsh.nel_per_level(ilev);
  if nel_lev == 0, continue; end
  u_lev = hspace.Csub{ilev} * u(1:ndofs);
  sp_tp = hspace.space_of_level(ilev);
  knx = sp_tp.knots{1};  kny = sp_tp.knots{2};
  dgx = sp_tp.degree(1); dgy = sp_tp.degree(2);
  coef = level_coefficients (hmsh, hspace, ilev, u_lev, sp_tp.ndof);
  n1 = numel (knx) - dgx - 1;  n2 = numel (kny) - dgy - 1;
  Cm = reshape (coef, n1, n2);
  [hxe, hye, cxe, cye] = elem_geom_of_level (hmsh, hspace, ilev, hmsh.active{ilev}(:));
  for iel = 1:nel_lev
    k = off + iel;
    ns = subcell_count (max (hxe(iel), hye(iel)), pb, q);
    [xq, wxq] = subcell_rule (cxe(iel)-hxe(iel)/2, cxe(iel)+hxe(iel)/2, ns, q.lay_nq);
    [yq, wyq] = subcell_rule (cye(iel)-hye(iel)/2, cye(iel)+hye(iel)/2, ns, q.lay_nq);
    [Bx, cix] = span_ders (knx, dgx, xq, 2);
    [By, ciy] = span_ders (kny, dgy, yq, 2);
    Cl = Cm(cix, ciy);
    V0  = Bx{1}*Cl*By{1}.';   Vx  = Bx{2}*Cl*By{1}.';
    Vy  = Bx{1}*Cl*By{2}.';   Vxx = Bx{3}*Cl*By{1}.';
    Vyy = Bx{1}*Cl*By{3}.';
    [XG, YG] = ndgrid (xq, yq);
    WG = wxq(:) * wyq(:).';
    RG = -pb.eps*(Vxx+Vyy) + pb.b(1)*Vx + pb.b(2)*Vy + pb.c*V0 - pb.f (XG,YG);
    r2(k) = max (sum (sum (WG .* RG.^2)), 0);
    e2(k) = max (sum (sum (WG .* (pb.u (XG,YG) - V0).^2)), 0);
  end
  off = off + nel_lev;
end
eta_ref = sqrt (sum (tauK.^2 .* r2));
err_ref = sqrt (sum (e2));
end


%% ======================================================================
%  CHARACTERISTIC SOLVE  (optional variant, PRM.COMPUTE_TR)
%  ======================================================================
function [etrK, dg] = transport_stream (LUT, pb, prm, ntot, NLforce)
%TRANSPORT_STREAM  eta_tr,K = ||etilde||_{L2(K)}, (b.grad+c)etilde=-R solved
%   along characteristics with local integrating factor, block-marched.
be = pb.ahat(:);  bn = pb.bnorm;
bp = [-be(2); be(1)];
cn = bp.' * [0 1 0 1; 0 0 1 1];
nmin = min (cn);  nmax = max (cn);
hmin = min ([LUT.hmin]);
if nargin < 5 || isempty (NLforce)
  NL = ceil (prm.tr_Ktr*(nmax-nmin)/max (hmin, eps));
  NL = max ([prm.tr_NLmin, 48, NL]);
  NL = min (NL, prm.tr_NLmax);
else
  NL = NLforce;
end
dn = (nmax - nmin) / NL;
nk = nmin + ((1:NL) - 0.5) * dn;

sxlo = ( 0 + nk*be(2))/be(1);  sxhi = ( 1 + nk*be(2))/be(1);
sylo = ( 0 - nk*be(1))/be(2);  syhi = ( 1 - nk*be(1))/be(2);
if be(1) < 0, t = sxlo; sxlo = sxhi; sxhi = t; end
if be(2) < 0, t = sylo; sylo = syhi; syhi = t; end
slo = max (sxlo, sylo);  shi = min (sxhi, syhi);
live = shi > slo + 1e-13;
slo = slo(live);  shi = shi(live);  nn = nk(live);  NLa = numel (nn);
if NLa == 0, etrK = zeros (ntot,1);  dg = empty_dg ();  return; end

etr2 = zeros (ntot,1);  area = 0;  nsegt = 0;  nstepmx = 0;  ndrop = 0;
blk = max (16, round (prm.tr_block));
for a0 = 1:blk:NLa
  z = a0 : min (a0+blk-1, NLa);
  [e2b, ab, nsb, nstb, ndb] = transport_block (LUT, pb, prm, ntot, ...
                                 nn(z), slo(z), shi(z), dn, be, bn);
  etr2 = etr2 + e2b;  area = area + ab;  nsegt = nsegt + nsb;
  nstepmx = max (nstepmx, nstb);  ndrop = ndrop + ndb;
end

etrK = sqrt (max (etr2, 0));
dg.NL = NL;  dg.nseg = nsegt;  dg.area = area;  dg.ndrop = ndrop;
dg.ppe = nsegt*prm.tr_nq/max (ntot,1);  dg.hmin = hmin;  dg.nstep = nstepmx;
if abs (dg.area - 1) > 2e-3 || ndrop > 0
  warning ('iga:area', ['swept area %.6f (should be 1) with %d unlocated ' ...
    'probes: eta_tr is NOT trustworthy on this mesh.'], dg.area, ndrop);
end
end

function [etr2, area, nseg, nstep, ndrop] = transport_block (LUT, pb, prm, ...
                                       ntot, nn, slo, shi, dn, be, bn)
cc = pb.c;
[xg, wg] = gauss_ref (prm.tr_nq);
Mc = cumgauss (prm.tr_nq);
dsfine = max (prm.tr_dsfac * pb.delta, 1e-12);
dsband = prm.tr_band * pb.delta;
xc1 = pb.xc(1);  xc2 = pb.xc(2);  r0 = pb.r0;
NLb = numel (nn);

% --- PASS 1 : lockstep march ----------------------------------------
snudge = 1e-11;  ndrop = 0;
s = slo(:).';  alive = shi > slo + 2*snudge;  step = 0;
maxstep = 200000;  CH = 256;
SA = zeros (CH, NLb);  DS = SA;  EL = SA;
while any (alive) && step < maxstep
  step = step + 1;
  if step > size (SA,1)
    SA(size (SA,1)+CH, NLb) = 0;
    DS(size (DS,1)+CH, NLb) = 0;
    EL(size (EL,1)+CH, NLb) = 0;
  end
  id  = find (alive);
  ssv = s(id).';  nnv = nn(id).';  shv = shi(id).';
  sm  = min (ssv + snudge, shv);
  X = sm*be(1) - nnv*be(2);
  Y = sm*be(2) + nnv*be(1);
  [gl, xlo, xhi, ylo, yhi] = locate_hier_batch (LUT, X, Y);
  dcirc = abs (hypot (X - xc1, Y - xc2) - r0);
  dsloc = min (prm.tr_dscoarse, dsfine + max (0, dcirc - dsband));
  if be(1) > 0,      sx = sm + (xhi - X)/be(1);
  elseif be(1) < 0,  sx = sm + (xlo - X)/be(1);
  else,              sx = inf (size (sm));
  end
  if be(2) > 0,      sy = sm + (yhi - Y)/be(2);
  elseif be(2) < 0,  sy = sm + (ylo - Y)/be(2);
  else,              sy = inf (size (sm));
  end
  s2 = min ([sx, sy, ssv + dsloc, shv], [], 2);
  s2 = max (s2, sm);
  ok = gl > 0;
  if any (~ok)                       % re-probe rather than drop a segment
    zb = find (~ok);
    Xb = X(zb) + 1e3*snudge*be(1);  Yb = Y(zb) + 1e3*snudge*be(2);
    [g2, xl2, xh2, yl2, yh2] = locate_hier_batch (LUT, Xb, Yb);
    gl(zb) = g2;
    if be(1) > 0,     sxb = sm(zb) + (xh2 - X(zb))/be(1);
    elseif be(1) < 0, sxb = sm(zb) + (xl2 - X(zb))/be(1);
    else,             sxb = inf (size (zb));
    end
    if be(2) > 0,     syb = sm(zb) + (yh2 - Y(zb))/be(2);
    elseif be(2) < 0, syb = sm(zb) + (yl2 - Y(zb))/be(2);
    else,             syb = inf (size (zb));
    end
    s2(zb) = max (min ([sxb, syb, ssv(zb) + dsloc(zb), shv(zb)], [], 2), sm(zb));
    ok = gl > 0;
    ndrop = ndrop + sum (~ok);
  end
  if any (ok)
    j = id(ok);
    SA(step, j) = ssv(ok).';
    DS(step, j) = (s2(ok) - ssv(ok)).';
    EL(step, j) = gl(ok).';
    s(j) = s2(ok).';
  end
  if any (~ok), s(id(~ok)) = s(id(~ok)) + snudge; end
  alive = s < shi - 2*snudge;
end
if step >= maxstep
  warning ('iga:march','streamline march hit the step cap; eta_tr is not trustworthy.');
end
nstep = step;
SA = SA(1:nstep,:);  DS = DS(1:nstep,:);  EL = EL(1:nstep,:);

% --- PASS 2 : batch residual evaluation ------------------------------
lin = find (DS > 0);   lin = lin(:);
[~, ci] = ind2sub (size (DS), lin);
ci = ci(:);
sA = SA(lin);  sA = sA(:);
dS = DS(lin);  dS = dS(:);
nseg = numel (lin);
area = sum (dS) * dn;
if nseg == 0, etr2 = zeros (ntot,1); return; end
Snod = sA(:, ones (1, prm.tr_nq)) + 0.5*dS*(1 + xg);
nnc  = nn(:);
Nnod = nnc(ci) * ones (1, prm.tr_nq);
Rn = hier_residual_scattered (LUT, pb, ...
       reshape (Snod*be(1) - Nnod*be(2), [], 1), ...
       reshape (Snod*be(2) + Nnod*be(1), [], 1));
Rn = reshape (Rn, nseg, prm.tr_nq);
clear Snod Nnod;
rowid = zeros (nstep, NLb);
rowid(lin) = 1:nseg;

% --- PASS 3 : sweep, vectorised across the tracks of this block -------
etilde = zeros (1, NLb);
accI = zeros (nseg,1);  accV = zeros (nseg,1);  ptr = 0;
for k = 1:nstep
  m = find (DS(k,:) > 0);
  if isempty (m), continue; end
  dm = DS(k,m).';  em = etilde(m).';  elm = EL(k,m).';
  Rk = Rn(rowid(k,m), :);
  tt = 0.5*dm*(1 + xg);
  psi = exp (cc*tt/bn) .* Rk;
  Gh  = 0.5*(dm*ones (1,prm.tr_nq)) .* (psi * Mc.');
  ev  = exp (-cc*tt/bn) .* (em*ones (1,prm.tr_nq) - Gh/bn);
  nm  = numel (m);
  accI(ptr+(1:nm)) = elm;
  accV(ptr+(1:nm)) = 0.5*dm .* (ev.^2 * wg(:)) * dn;
  ptr = ptr + nm;
  Ge = 0.5*dm .* (psi * wg(:));
  etilde(m) = (exp (-cc*dm/bn) .* (em - Ge/bn)).';
end
etr2 = accumarray (accI(1:ptr), accV(1:ptr), [ntot 1]);

if prm.bl_correction
  lastk = zeros (1, NLb);
  for k = 1:nstep, mm = DS(k,:) > 0;  lastk(mm) = k; end
  gk = find (lastk > 0);
  if ~isempty (gk)
    eidx = EL (sub2ind (size (EL), lastk(gk), gk));
    wbl  = (pb.eps/(2*bn)) * (etilde(gk).^2) * dn;
    etr2 = etr2 + accumarray (eidx(:), wbl(:), [ntot 1]);
  end
end
end

function s = gradind_scale_factor (hx, hy, prm)
if strcmpi (prm.gradind_scale, 'h'), s = sqrt (hx(:) .* hy(:));
else,                                s = ones (numel (hx), 1); end
end

%% ======================================================================
%  SELF TEST
%  ======================================================================
function out = run_selftest (prm)
fprintf ('\n  SELF TEST\n  %s\n', repmat ('=',1,70));
ok = true;

for nq = [3 5 6 8]
  M = cumgauss (nq);  [xg,~] = gauss_ref (nq);
  e = max (abs (sum (M,2).' - (xg+1)));
  fprintf ('  cumgauss nq=%d   rowsum identity  %.2e   %s\n', nq, e, pf (e<1e-12));
  ok = ok && e < 1e-12;
end

[CA,CD] = vms_constants (setfield (prm, 'norm_type', 'hs'));   %#ok<SFLD>
e1 = abs (CA - 1/sqrt (2)) + abs (CD - 1/(3*sqrt (10)));
[CA,CD] = vms_constants (setfield (prm, 'norm_type', 'operator'));  %#ok<SFLD>
e2 = abs (CA - 2/pi) + abs (CD - 1/pi^2);
fprintf ('  vms_constants (paper 4.7/4.8)     %.2e   %s\n', e1+e2, pf (e1+e2<1e-14));
ok = ok && (e1+e2) < 1e-14;

pb = ex53_problem (prm);
[nu, ng] = exact_norms (pb);
fprintf ('  ||u_exact||_L2   = %.6f   (reference 0.363482)\n', nu);
fprintf ('  ||grad u||_L2    = %.6f   (reference 4.068946)\n', ng);
ok = ok && abs (nu-0.363482) < 1e-5 && abs (ng-4.068946) < 1e-4;

t = linspace (0,1,2001).';
v = max (abs ([pb.u(0*t,t); pb.u(0*t+1,t); pb.u(t,0*t); pb.u(t,0*t+1)]));
fprintf ('  max|u| on dOmega = %.2e                       %s\n', v, pf (v<1e-12));
ok = ok && v < 1e-12;

h = 1e-5;  rng ('default');
xs = 0.1+0.8*rand (200,1);  ys = 0.1+0.8*rand (200,1);
xs(1:20) = 0.5 + 0.25*cos (linspace (0,2*pi,20).');
ys(1:20) = 0.5 + 0.25*sin (linspace (0,2*pi,20).');
lap = (pb.u(xs+h,ys)+pb.u(xs-h,ys)+pb.u(xs,ys+h)+pb.u(xs,ys-h)-4*pb.u(xs,ys))/h^2;
ux  = (pb.u(xs+h,ys)-pb.u(xs-h,ys))/(2*h);
uy  = (pb.u(xs,ys+h)-pb.u(xs,ys-h))/(2*h);
fd  = -pb.eps*lap + pb.b(1)*ux + pb.b(2)*uy + pb.c*pb.u(xs,ys);
r   = max (abs (fd - pb.f(xs,ys)) ./ max (abs (pb.f(xs,ys)),1));
fprintf ('  f vs finite differences           %.2e   %s  (FD truncation)\n', r, pf (r<1e-3));

fprintf ('  %s\n', repmat ('-',1,70));
[rel, ar] = transport_identity_test (prm);
fprintf ('  transport identity  rel.err %.2e  swept area %.8f   %s\n', ...
         rel, ar, pf (rel < 5e-4 && abs (ar-1) < 1e-4));
ok = ok && rel < 5e-4 && abs (ar-1) < 1e-4;

fprintf ('  %s\n  SELF TEST %s\n\n', repmat ('=',1,70), tern (ok,'PASSED','FAILED'));
out = struct ('ok', ok);
end

function [rel, ar] = transport_identity_test (prm)
%TRANSPORT_IDENTITY_TEST  eps=0, f=0 => etilde=-u_h exactly; checks eta_tr vs ||u_h||.
p = 4; reg = 3; nsub = 8;
geometry = geo_load (nrb4surf ([0 0],[1 0],[0 1],[1 1]));
[knots, zeta] = kntrefine (geometry.nurbs.knots, (nsub-1)*ones(1,2), ...
                           p*ones(1,2), reg*ones(1,2));
rule = msh_gauss_nodes ((p+4)*ones(1,2));
[qn,qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, geometry);
space = sp_bspline (knots, p*ones(1,2), msh);
hmsh   = hierarchical_mesh (msh,[2 2]);
hspace = hierarchical_space (hmsh, space, 'simplified', true, reg*ones(1,2));
rng ('default');
nd = hspace.ndof;  u = randn (nd,1);
n1 = round (sqrt (nd));
U = reshape (u, n1, n1);  U(1,:) = 0;  U(:,1) = 0;   % zero on the inflow sides
u = U(:);
pbz = ex53_problem (prm);
pbz.eps = 0;  pbz.f = @(x,y) zeros (size (x));
q = prm;  q.tr_NLmin = 356;  q.tr_nq = 8;  q.tr_dsfac = 0.5;
q.tr_dscoarse = 0.05;  q.tr_block = 128;  q.bl_correction = false;
LUT = build_level_lut (hmsh, hspace);
LUT = attach_coefficients (LUT, hmsh, hspace, u);
[etrK, dg] = transport_stream (LUT, pbz, q, hmsh.nel, []);
eta = sqrt (sum (etrK.^2));
ex = 0;
for ilev = 1:hmsh.nlevels
  nel = hmsh.nel_per_level(ilev);
  if nel == 0, continue; end
  [msh_c, sp_c] = eval_elems (hmsh, hspace, ilev, hmsh.active{ilev}(:));
  ul = hspace.Csub{ilev} * u(1:hspace.ndof_per_level(ilev));
  for iel = 1:nel
    w = msh_c.quad_weights(:,iel) .* msh_c.jacdet(:,iel);
    N = reshape (sp_c.shape_functions(:,:,iel), msh_c.nqn, sp_c.nsh_max);
    vh = N*ul(sp_c.connectivity(:,iel));
    ex = ex + sum (w .* vh.^2);
  end
end
ex = sqrt (ex);
rel = abs (eta/ex - 1);  ar = dg.area;
end

function s = pf (t)
if t, s = 'ok'; else, s = 'FAIL'; end
end
function s = tern (t,a,b)
if t, s = a; else, s = b; end
end

function [nu, ng] = exact_norms (pb)
[xg, wg] = gauss_ref (10);
e = linspace (0,1,300);
xs = zeros (1,(numel(e)-1)*10);  ws = xs;
for i = 1:numel (e)-1
  z = (i-1)*10 + (1:10);
  xs(z) = 0.5*(e(i)+e(i+1)) + 0.5*(e(i+1)-e(i))*xg;
  ws(z) = 0.5*(e(i+1)-e(i))*wg;
end
nu = 0;  ng = 0;
for j = 1:numel (xs)
  [uu, ux, uy] = pb.ug (xs.', xs(j)*ones (numel (xs),1));
  nu = nu + ws(j)*sum (ws(:).*uu.^2);
  ng = ng + ws(j)*sum (ws(:).*(ux.^2+uy.^2));
end
nu = sqrt (nu);  ng = sqrt (ng);
end


%% ======================================================================
%  DAMKOHLER SWEEP
%  ======================================================================
function out = run_damkohler (prm)
fprintf ('\n  DAMKOHLER SWEEP   (|b| = %.4f, diam(Omega) = %.4f)\n', ...
         norm (prm.b), sqrt (2));
fprintf ('  %s\n', repmat ('=',1,76));
fprintf ('  %6s %10s %8s | %4s %8s | %11s %11s | %8s\n', ...
         'c','L_decay','L/diam','p','ndof','||e||','eta','I_eff');
fprintf ('  %s\n', repmat ('-',1,76));
q = prm;  q.verbose = false;  q.damkohler = [];  q.refinement = 'uniform';
q.plots = struct ('error',false,'ieff',false,'maps',false,'mesh',false,'solution',false);
rows = {};
for c = prm.damkohler(:).'
  q.c = c;
  for p = prm.degrees(:).'
    q.degrees = p;
    R = run_sweep (q, 'uniform');
    if isempty (R) || ~R(1).ok, continue; end
    k = numel (R(1).ndof);
    Ld = norm (prm.b)/c;
    fprintf ('  %6.1f %10.4f %8.3f | %4d %8d | %11.4e %11.4e | %8.3f\n', ...
             c, Ld, Ld/sqrt (2), p, R(1).ndof(k), R(1).err(k), ...
             R(1).eta(k), R(1).ieff(k));
    rows{end+1} = struct ('c',c,'p',p,'R',R(1));  %#ok<AGROW>
  end
end
fprintf ('  %s\n', repmat ('-',1,76));
out = struct ('rows', {rows}, 'prm', prm);
end


%% ======================================================================
%  OPTION HANDLING
%  ======================================================================
function prm = validate_options (prm)
if ischar (prm.drivers), prm.drivers = {prm.drivers}; end
for i = 1:numel (prm.drivers)
  if ~any (strcmpi (prm.drivers{i}, {'vms','gradind'}))
    error ('iga:driver','Unknown driver ''%s''.', prm.drivers{i});
  end
end
if ischar (prm.fig_fmt), prm.fig_fmt = {prm.fig_fmt}; end
if ~any (strcmpi (prm.gradind_scale, {'none','h'}))
  error ('iga:gradscale','gradind_scale must be ''none'' or ''h''.');
end
if ~any (strcmpi (prm.mark_on, {'transport','local'}))
  error ('iga:markon','mark_on must be ''transport'' or ''local''.');
end
if ~any (strcmpi (prm.space_type, {'standard','simplified'}))
  error ('iga:spacetype','space_type must be ''standard'' or ''simplified''.');
end
if strcmpi (prm.mark_on,'transport') && ~prm.compute_tr
  error ('iga:markon','mark_on = ''transport'' needs compute_tr = true.');
end
if strcmpi (prm.est_primary,'transport') && ~prm.compute_tr
  error ('iga:prim','est_primary = ''transport'' needs compute_tr = true.');
end
if prm.gamma_safe ~= 1
  warning ('iga:gamma', ['gamma_safe = %.3f.  The theory licenses no such ' ...
    'constant; it moves the offset of I_eff and never its spread.'], prm.gamma_safe);
end
end

function drivers = active_drivers (prm)
if strcmpi (prm.refinement, 'uniform'), drivers = {'uniform'};
else,                                   drivers = prm.drivers; end
end

function prm = parse_options (prm, args)
if numel (args) == 1 && isstruct (args{1})
  f = fieldnames (args{1});
  for i = 1:numel (f), prm.(f{i}) = args{1}.(f{i}); end
  return;
end
if mod (numel (args), 2) ~= 0
  error ('iga:args','Options must come in name/value pairs.');
end
dead = {'combine','tr_grid','estimator','estimators','indicator','tau_rule'};
for i = 1:2:numel (args)
  name = args{i};
  if ~ischar (name), error ('iga:args','Option names must be strings.'); end
  if any (strcmpi (name, dead))
    error ('iga:args','''%s'' no longer exists.  Use ''est_primary''.', name);
  end
  if isfield (prm, name),            prm.(name) = args{i+1};
  elseif isfield (prm.plots, name),  prm.plots.(name) = args{i+1};
  else,  error ('iga:args','Unknown option ''%s''.', name);
  end
end
end


%% ======================================================================
%  SWEEP / DRIVER
%  ======================================================================
function [R, pb, niter] = run_sweep (prm, driver)
if nargin < 2, driver = 'vms'; end
SWEEP = [1 0; 2 1; 3 2; 4 3; 5 4];
prm.degrees = prm.degrees(:).';
if any (prm.degrees < 1 | prm.degrees > 5)
  error ('iga:degrees','prm.degrees must be a subset of 1:5.');
end
if isempty (prm.niter)
  if strcmpi (prm.refinement,'uniform'), niter = prm.niter_uniform;
  else,                                  niter = prm.niter_adaptive; end
else
  niter = prm.niter;
end
pb = ex53_problem (prm);
check_boundary_data (pb);
[~, pb.gradref] = exact_norms (pb);
R = struct ('degree',{},'regularity',{},'driver',{},'ndof',{},'err',{}, ...
            'eta',{},'ieff',{},'eta_loc',{},'eta_tr',{},'ieff_loc',{}, ...
            'ieff_tr',{},'eta_grad',{},'ieff_grad',{},'eta_cert',{}, ...
            'ieff_cert',{},'rcond',{},'area',{},'ppe',{},'NL',{}, ...
            'efloor',{},'trust',{},'tsec',{}, ...
            'nel',{},'hmin',{},'err_proj',{},'rstab',{}, ...
            'fel',{},'fet',{},'eta_ref',{},'err_ref',{},'raw',{}, ...
            'rho_v',{},'rho_g',{},'cap_v',{},'cap_g',{},'satg',{},'aspect',{}, ...
            'tr_rich',{},'loc',{},'mq_loc',{},'mq_grad',{},'sol',{},'ok',{});
for k = 1:numel (prm.degrees)
  id = prm.degrees(k);
  p = SWEEP(id,1);  reg = SWEEP(id,2);
  R(k).degree = p;  R(k).regularity = reg;  R(k).driver = driver;  R(k).ok = false;
  R(k).mq_loc = empty_mq (); R(k).mq_grad = empty_mq ();
  try
    S = run_one_degree (p, reg, niter, pb, prm, driver);
    f = fieldnames (S);
    for j = 1:numel (f), R(k).(f{j}) = S.(f{j}); end
    R(k).ok = true;
  catch err
    fprintf (2, '  p = %d, driver ''%s'' FAILED: %s\n', p, driver, err.message);
    if prm.verbose && ~isempty (err.stack)
      for q = 1:min (3, numel (err.stack))
        fprintf (2, '     at %s line %d\n', err.stack(q).name, err.stack(q).line);
      end
    end
  end
end
end

function r = fit_rate (nd, er)
g = isfinite (er) & er > 0 & isfinite (nd) & nd > 0;
nd = nd(g);  er = er(g);
if numel (nd) < 2, r = [NaN NaN]; return; end
z = max (1, numel (nd)-4) : numel (nd);
r = polyfit (log (nd(z)), log (er(z)), 1);
end
function v = sub_first (r)
v = r(1);
end
function m = empty_mq ()
m = struct ('rho_rank',NaN,'capture',NaN,'frac_marked',NaN,'frac_ideal',NaN);
end

function S = run_one_degree (p, reg, niter, pb, prm, driver)
q = prm;  q.degree = p;  q.regularity = reg;
if isempty (q.nquad), q.nquad = p + 2; end
geometry = geo_load (nrb4surf ([0 0],[1 0],[0 1],[1 1]));
[knots, zeta] = kntrefine (geometry.nurbs.knots, (q.nsub0-1)*ones (1,2), ...
                           p*ones (1,2), reg*ones (1,2));
rule     = msh_gauss_nodes (q.nquad*ones (1,2));
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);
space    = sp_bspline (knots, p*ones (1,2), msh);
hmsh   = hierarchical_mesh (msh, [2 2]);
hspace = hierarchical_space (hmsh, space, q.space_type, q.truncated, ...
                             reg*ones (1,2));

f = {'ndof','err','eta','ieff','eta_loc','eta_tr','ieff_loc','ieff_tr', ...
     'eta_grad','ieff_grad','eta_cert','ieff_cert','rcond','area','ppe', ...
     'NL','tr_rich','efloor','trust','tsec','nel','hmin','err_proj','rstab', ...
     'fel','fet','eta_ref','err_ref','rho_v','rho_g','cap_v','cap_g', ...
     'satg','aspect'};
S.raw = {};
for i = 1:numel (f), S.(f{i}) = []; end
E = [];  u = [];
for it = 1:niter
  t0 = tic;
  [u, rc, relcorr] = solve_stabilised (hmsh, hspace, pb, q);
  E = estimate_and_error (hmsh, hspace, u, pb, q);
  E.eta_mark = mark_vector (E, driver, q);
  S.ndof(end+1)  = hspace.ndof;
  S.err(end+1)   = E.err;
  S.eta(end+1)   = E.eta;      S.ieff(end+1)   = E.eta / max (E.err, realmin);
  S.eta_loc(end+1)= E.eta_loc; S.ieff_loc(end+1)= E.eta_loc/max (E.err,realmin);
  S.eta_tr(end+1) = E.eta_tr;  S.ieff_tr(end+1) = E.eta_tr /max (E.err,realmin);
  S.eta_grad(end+1)= E.eta_grad; S.ieff_grad(end+1)= E.eta_grad/max (E.err,realmin);
  S.eta_cert(end+1)= E.eta_cert; S.ieff_cert(end+1)= E.eta_cert/max (E.err,realmin);
  S.rcond(end+1) = rc;
  S.efloor(end+1) = relcorr * E.unorm;
  S.trust(end+1)  = E.err > q.err_floor_factor * S.efloor(end);
  S.area(end+1)  = E.dg.area;  S.ppe(end+1) = E.dg.ppe;  S.NL(end+1) = E.dg.NL;
  S.tr_rich(end+1) = E.tr_rich;
  S.nel(end+1)  = hmsh.nel;
  S.hmin(end+1) = min (min (E.hx), min (E.hy));
  S.raw{end+1}  = struct ('hx', E.hx, 'hy', E.hy, 'r', E.rK, 'e', E.err_K);
  S.fel(end+1)  = E.frac_err_lay;   S.fet(end+1) = E.frac_eta_lay;
  S.eta_ref(end+1) = E.eta_ref;     S.err_ref(end+1) = E.err_ref;
  % ---- scale-free indicator quality (kept in OUT, not printed) --------
  S.rho_v(end+1) = pearson_ (rank_ (E.eta_locK),  rank_ (E.err_K));
  S.rho_g(end+1) = pearson_ (rank_ (E.eta_gradK), rank_ (E.err_K));
  S.cap_v(end+1) = capture_ratio (E.eta_locK,  E.err_K, q.theta);
  S.cap_g(end+1) = capture_ratio (E.eta_gradK, E.err_K, q.theta);
  S.satg(end+1)  = E.eta_grad / max (pb.gradref, realmin);
  S.aspect(end+1)= max (max (E.hx), max (E.hy)) / ...
                   max (min (min (E.hx), min (E.hy)), realmin);
  if q.diag_proj
    S.err_proj(end+1) = l2_projection_error (hmsh, hspace, pb, q);
    S.rstab(end+1)    = E.err / max (S.err_proj(end), realmin);
  else
    S.err_proj(end+1) = NaN;  S.rstab(end+1) = NaN;
  end
  S.tsec(end+1) = toc (t0);
  if q.verbose            % report the estimator that is DRIVING this run
    if strcmpi (driver, 'gradind')
      ed = S.eta_grad(end);  ei = S.ieff_grad(end);
    else
      ed = S.eta(end);       ei = S.ieff(end);
    end
    fprintf (['    p=%d %-8s it %2d   ndof %6d  nel %6d   ||e|| %.4e  ' ...
              'eta %.4e  I_eff %9.4f   %.1fs\n'], p, driver, it, ...
             hspace.ndof, hmsh.nel, E.err, ed, ei, S.tsec(end));
  end
  if it == niter, break; end
  if q.stop_at_floor && ~S.trust(end)
    if q.verbose
      fprintf ('    [stop: ||e|| %.3e has reached the solver floor %.3e]\n', ...
               E.err, S.efloor(end));
    end
    break;
  end
  [hmsh2, hspace2, moved] = refine_once (hmsh, hspace, E, q);
  if ~moved, break; end
  if hspace2.ndof > q.max_ndof, break; end
  hmsh = hmsh2;  hspace = hspace2;
end
S.loc     = local_effectivity (E, E.eta_K, pb, q);
S.mq_loc  = marker_quality (E.eta_locK,  E.err_K, q.theta);
S.mq_grad = marker_quality (E.eta_gradK, E.err_K, q.theta);
S.sol.hmsh = hmsh;  S.sol.hspace = hspace;  S.sol.u = u;  S.sol.geometry = geometry;
end

function [hmsh, hspace, moved] = refine_once (hmsh, hspace, E, q)
moved = false;
use_fun = strcmpi (q.mark_object,'functions') && ~strcmpi (q.refinement,'uniform');
if use_fun
  try
    mf = mark_functions (hmsh, hspace, E.eta_mark, q);
    if sum (cellfun (@numel, mf)) > 0
      adap.flag = 'functions';
      [hmsh, hspace] = adaptivity_refine (hmsh, hspace, mf, adap);
      moved = true;  return;
    end
  catch merr
    warning ('iga:markfun','function marking unavailable (%s); using elements.', merr.message);
  end
end
marked = mark_elements (hmsh, E.eta_mark, E.lev_of, E.pos_of, q);
if q.adm_class >= 2
  marked = enforce_admissibility (hmsh, hspace, marked, q);
end
if sum (cellfun (@numel, marked)) == 0, return; end
adap.flag = 'elements';
[hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adap);
moved = true;
end

function v = mark_vector (E, driver, prm)
switch lower (driver)
  case 'gradind', v = E.eta_gradK;
  otherwise
    if strcmpi (prm.mark_on,'local'), v = E.eta_locK; else, v = E.eta_trK; end
end
end

%% ======================================================================
%  SOLVER  (chunked assembly, equilibration, iterative refinement)
%  ======================================================================
function [u, rc, relcorr] = solve_stabilised (hmsh, hspace, pb, prm)
II = {};  JJ = {};  VV = {};
rhs = zeros (hspace.ndof, 1);
ndofs = 0;
for ilev = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(ilev);
  nel_lev = hmsh.nel_per_level(ilev);
  if nel_lev == 0, continue; end
  act = hmsh.active{ilev}(:);
  nrows = size (hspace.Csub{ilev}, 1);
  Ic = {};  Jc = {};  Vc = {};
  Floc = zeros (nrows, 1);
  sp_tp = hspace.space_of_level(ilev);
  knx = sp_tp.knots{1};   kny = sp_tp.knots{2};
  dgx = sp_tp.degree(1);  dgy = sp_tp.degree(2);
  do_sub = strcmpi (prm.quad_mode, 'subcell');
  band   = prm.lay_bandw * pb.delta;
  nc = chunk_size (hspace, ilev, prm.nquad^2, prm);
  for a0 = 1:nc:nel_lev
    sel = a0 : min (a0+nc-1, nel_lev);
    [msh_c, sp_c] = eval_elems (hmsh, hspace, ilev, act(sel));
    [hxe, hye, cxe, cye] = elem_geom_of_level (hmsh, hspace, ilev, act(sel));
    nsh = sp_c.nsh_max;  nqn = msh_c.nqn;  nel = numel (sel);
    I = zeros (nsh*nsh*nel,1);  J = I;  V = I;  ptr = 0;
    for iel = 1:nel
      w = msh_c.quad_weights(:,iel) .* msh_c.jacdet(:,iel);
      x = reshape (msh_c.geo_map(1,:,iel), nqn, 1);
      y = reshape (msh_c.geo_map(2,:,iel), nqn, 1);
      [N, Gx, Gy, Lp] = shape_blocks (sp_c, iel, nqn, nsh);
      bG  = pb.b(1)*Gx + pb.b(2)*Gy;
      Rop = -pb.eps*Lp + bG + pb.c*N;
      if prm.stabilise
        hf   = elem_length (hxe(iel), hye(iel), pb.ahat, prm.h_stab);
        alph = pb.bnorm * hf / (2 * pb.eps);
        tau  = hf / (2 * pb.bnorm) * xi_coth (alph);       % paper (2.26)
        switch lower (prm.test_op)
          case 'sgs',  P = bG + pb.eps*Lp - pb.c*N;        % -L^* v  (MS)
          case 'supg', P = bG;
          case 'gls',  P = Rop;
          otherwise, error ('iga:testop','Unknown test_op ''%s''.', prm.test_op);
        end
      else
        tau = 0;  P = zeros (nqn, nsh);
      end
      fq = pb.f (x, y);
      W  = w(:, ones (1, nsh));
      Ke = N.'*(W.*bG) + pb.eps*(Gx.'*(W.*Gx) + Gy.'*(W.*Gy)) ...
           + pb.c*(N.'*(W.*N)) + tau*(P.'*(W.*Rop));
      if do_sub && circle_crosses (cxe(iel), cye(iel), hxe(iel), hye(iel), pb, band)
        Fe = elem_load_subcell (knx, kny, dgx, dgy, ...
               cxe(iel)-hxe(iel)/2, cxe(iel)+hxe(iel)/2, ...
               cye(iel)-hye(iel)/2, cye(iel)+hye(iel)/2, ...
               max (hxe(iel), hye(iel)), tau, pb, prm);
      else
        Fe = N.'*(w.*fq) + tau*(P.'*(w.*fq));
      end
      conn = sp_c.connectivity(:,iel);
      idx = ptr + (1:nsh*nsh);
      I(idx) = repmat (conn, nsh, 1);
      J(idx) = reshape (repmat (conn.', nsh, 1), [], 1);
      V(idx) = Ke(:);
      ptr = ptr + nsh*nsh;
      Floc(conn) = Floc(conn) + Fe;
    end
    Ic{end+1} = I(1:ptr);  Jc{end+1} = J(1:ptr);  Vc{end+1} = V(1:ptr); %#ok<AGROW>
    clear msh_c sp_c I J V;
  end
  A_TP = sparse (cat (1,Ic{:}), cat (1,Jc{:}), cat (1,Vc{:}), nrows, nrows);
  clear Ic Jc Vc;
  C = hspace.Csub{ilev};
  B = C.' * A_TP * C;
  [bi, bj, bv] = find (B);
  II{end+1} = bi;  JJ{end+1} = bj;  VV{end+1} = bv;   %#ok<AGROW>
  rhs(1:ndofs) = rhs(1:ndofs) + C.' * Floc;
end
A = sparse (cat (1, II{:}), cat (1, JJ{:}), cat (1, VV{:}), ...
            hspace.ndof, hspace.ndof);
clear II JJ VV;
drchlt = dirichlet_dofs (hmsh, hspace);
free   = setdiff (1:hspace.ndof, drchlt);
Af = A(free,free);  bf = rhs(free);
clear A;
% --- symmetric diagonal equilibration (levels differ in scale by orders of magnitude)
dv = full (sqrt (abs (diag (Af))));
dv(~isfinite (dv) | dv <= 0) = 1;
nf = numel (dv);
D  = spdiags (1./dv, 0, nf, nf);
As = D*Af*D;   bs = D*bf;
[L_,U_,P_,Q_,Rr] = lu (As);
xs = Q_ * (U_ \ (L_ \ (P_ * (Rr \ bs))));
% --- iterative refinement; correction size defines the accuracy floor
relcorr = 0;
for kir = 1:max (1, prm.iter_refine)
  rs = bs - As*xs;
  dx = Q_ * (U_ \ (L_ \ (P_ * (Rr \ rs))));
  xs = xs + dx;
  relcorr = norm (dx) / max (norm (xs), realmin);
end
u = zeros (hspace.ndof, 1);
u(free) = D*xs;
rc = NaN;
if prm.do_condest && nf <= prm.condest_max
  try, rc = 1/condest (As); catch, rc = NaN; end
end
end


function [eproj, uP] = l2_projection_error (hmsh, hspace, pb, prm)
%L2_PROJECTION_ERROR  ||u - Pi_h u||_L2 onto same THB space (off by default).
II = {};  JJ = {};  VV = {};
rhs = zeros (hspace.ndof, 1);
ndofs = 0;
band = prm.lay_bandw * pb.delta;
do_sub = strcmpi (prm.quad_mode, 'subcell');
for ilev = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(ilev);
  nel_lev = hmsh.nel_per_level(ilev);
  if nel_lev == 0, continue; end
  act = hmsh.active{ilev}(:);
  nrows = size (hspace.Csub{ilev}, 1);
  sp_tp = hspace.space_of_level(ilev);
  knx = sp_tp.knots{1};   kny = sp_tp.knots{2};
  dgx = sp_tp.degree(1);  dgy = sp_tp.degree(2);
  Ic = {};  Jc = {};  Vc = {};
  Floc = zeros (nrows, 1);
  nc = chunk_size (hspace, ilev, prm.nquad^2, prm);
  for a0 = 1:nc:nel_lev
    sel = a0 : min (a0+nc-1, nel_lev);
    [msh_c, sp_c] = eval_elems (hmsh, hspace, ilev, act(sel));
    [hxe, hye, cxe, cye] = elem_geom_of_level (hmsh, hspace, ilev, act(sel));
    nsh = sp_c.nsh_max;  nqn = msh_c.nqn;  nel = numel (sel);
    I = zeros (nsh*nsh*nel,1);  J = I;  V = I;  ptr = 0;
    for iel = 1:nel
      w = msh_c.quad_weights(:,iel) .* msh_c.jacdet(:,iel);
      x = reshape (msh_c.geo_map(1,:,iel), nqn, 1);
      y = reshape (msh_c.geo_map(2,:,iel), nqn, 1);
      N = reshape (sp_c.shape_functions(:,:,iel), nqn, nsh);
      W = w(:, ones (1, nsh));
      Me = N.'*(W.*N);
      if do_sub && circle_crosses (cxe(iel), cye(iel), hxe(iel), hye(iel), pb, band)
        T = elem_moment_subcell (knx, kny, dgx, dgy, ...
              cxe(iel)-hxe(iel)/2, cxe(iel)+hxe(iel)/2, ...
              cye(iel)-hye(iel)/2, cye(iel)+hye(iel)/2, ...
              max (hxe(iel), hye(iel)), pb.u, pb, prm);
        Fe = T(:);
      else
        Fe = N.'*(w.*pb.u (x,y));
      end
      conn = sp_c.connectivity(:,iel);
      idx = ptr + (1:nsh*nsh);
      I(idx) = repmat (conn, nsh, 1);
      J(idx) = reshape (repmat (conn.', nsh, 1), [], 1);
      V(idx) = Me(:);
      ptr = ptr + nsh*nsh;
      Floc(conn) = Floc(conn) + Fe;
    end
    Ic{end+1} = I(1:ptr);  Jc{end+1} = J(1:ptr);  Vc{end+1} = V(1:ptr); %#ok<AGROW>
    clear msh_c sp_c I J V;
  end
  M_TP = sparse (cat (1,Ic{:}), cat (1,Jc{:}), cat (1,Vc{:}), nrows, nrows);
  C = hspace.Csub{ilev};
  B = C.' * M_TP * C;
  [bi, bj, bv] = find (B);
  II{end+1} = bi;  JJ{end+1} = bj;  VV{end+1} = bv;   %#ok<AGROW>
  rhs(1:ndofs) = rhs(1:ndofs) + C.' * Floc;
end
M = sparse (cat (1, II{:}), cat (1, JJ{:}), cat (1, VV{:}), ...
            hspace.ndof, hspace.ndof);
drchlt = dirichlet_dofs (hmsh, hspace);
free   = setdiff (1:hspace.ndof, drchlt);
uP = zeros (hspace.ndof, 1);
uP(free) = M(free,free) \ rhs(free);
eproj = l2_err_of (hmsh, hspace, uP, pb, prm);
end

function T = elem_moment_subcell (knx, kny, px, py, x0, x1, y0, y1, hmax, ...
                                  fun, pb, prm)
%ELEM_MOMENT_SUBCELL  (fun, N_A)_K via layer sub-cell rule.
ns = subcell_count (hmax, pb, prm);
[xq, wxq] = subcell_rule (x0, x1, ns, prm.lay_nq);
[yq, wyq] = subcell_rule (y0, y1, ns, prm.lay_nq);
Bx = span_ders (knx, px, xq, 0);
By = span_ders (kny, py, yq, 0);
[XG, YG] = ndgrid (xq, yq);
Wf = (wxq(:) * wyq(:).') .* fun (XG, YG);
T = Bx{1}.' * Wf * By{1};
end

function er = l2_err_of (hmsh, hspace, u, pb, prm)
%L2_ERR_OF  ||u_exact - u_h||_L2 with the same quadrature as ESTIMATE_AND_ERROR.
e2 = 0;  ndofs = 0;
sub  = strcmpi (prm.quad_mode, 'subcell');
band = prm.lay_bandw * pb.delta;
for ilev = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(ilev);
  nel_lev = hmsh.nel_per_level(ilev);
  if nel_lev == 0, continue; end
  u_lev = hspace.Csub{ilev} * u(1:ndofs);
  sp_tp = hspace.space_of_level(ilev);
  knx = sp_tp.knots{1};  kny = sp_tp.knots{2};
  dgx = sp_tp.degree(1); dgy = sp_tp.degree(2);
  if sub
    coef = level_coefficients (hmsh, hspace, ilev, u_lev, sp_tp.ndof);
    n1 = numel (knx) - dgx - 1;  n2 = numel (kny) - dgy - 1;
    Cm = reshape (coef, n1, n2);
  end
  act = hmsh.active{ilev}(:);
  nc  = chunk_size (hspace, ilev, prm.nquad^2, prm);
  for a0 = 1:nc:nel_lev
    sel = a0 : min (a0+nc-1, nel_lev);
    [msh_c, sp_c] = eval_elems (hmsh, hspace, ilev, act(sel));
    [hxe, hye, cxe, cye] = elem_geom_of_level (hmsh, hspace, ilev, act(sel));
    nsh = sp_c.nsh_max;  nqn = msh_c.nqn;
    for iel = 1:numel (sel)
      w = msh_c.quad_weights(:,iel) .* msh_c.jacdet(:,iel);
      x = reshape (msh_c.geo_map(1,:,iel), nqn, 1);
      y = reshape (msh_c.geo_map(2,:,iel), nqn, 1);
      N = reshape (sp_c.shape_functions(:,:,iel), nqn, nsh);
      uh = N*u_lev(sp_c.connectivity(:,iel));
      if sub && circle_crosses (cxe(iel), cye(iel), hxe(iel), hye(iel), pb, band)
        ns = subcell_count (max (hxe(iel), hye(iel)), pb, prm);
        [xq, wxq] = subcell_rule (cxe(iel)-hxe(iel)/2, cxe(iel)+hxe(iel)/2, ns, prm.lay_nq);
        [yq, wyq] = subcell_rule (cye(iel)-hye(iel)/2, cye(iel)+hye(iel)/2, ns, prm.lay_nq);
        [Bx, cix] = span_ders (knx, dgx, xq, 0);
        [By, ciy] = span_ders (kny, dgy, yq, 0);
        V0 = Bx{1}*Cm(cix,ciy)*By{1}.';
        [XG, YG] = ndgrid (xq, yq);
        WG = wxq(:) * wyq(:).';
        e2 = e2 + max (sum (sum (WG .* (pb.u (XG,YG) - V0).^2)), 0);
      else
        e2 = e2 + max (sum (w .* (pb.u (x,y) - uh).^2), 0);
      end
    end
    clear msh_c sp_c;
  end
end
er = sqrt (max (e2, 0));
end

function Fe = elem_load_subcell (knx, kny, px, py, x0, x1, y0, y1, hmax, ...
                                 tau, pb, prm)
%ELEM_LOAD_SUBCELL  Fe = (N,f)_K + tau*(PN,f)_K via sub-cell rule on the layer.
ns = subcell_count (hmax, pb, prm);
[xq, wxq] = subcell_rule (x0, x1, ns, prm.lay_nq);
[yq, wyq] = subcell_rule (y0, y1, ns, prm.lay_nq);
Bx = span_ders (knx, px, xq, 2);
By = span_ders (kny, py, yq, 2);
[XG, YG] = ndgrid (xq, yq);
Wf = (wxq(:) * wyq(:).') .* pb.f (XG, YG);
T00 = Bx{1}.' * Wf * By{1};
T10 = Bx{2}.' * Wf * By{1};
T01 = Bx{1}.' * Wf * By{2};
T20 = Bx{3}.' * Wf * By{1};
T02 = Bx{1}.' * Wf * By{3};
Tlap = T20 + T02;
Tadv = pb.b(1)*T10 + pb.b(2)*T01;
switch lower (prm.test_op)
  case 'sgs',  TP = Tadv + pb.eps*Tlap - pb.c*T00;
  case 'supg', TP = Tadv;
  case 'gls',  TP = -pb.eps*Tlap + Tadv + pb.c*T00;
  otherwise, error ('iga:testop','Unknown test_op ''%s''.', prm.test_op);
end
Fe = T00 + tau*TP;
Fe = Fe(:);
end


%% ======================================================================
%  MARKING, ADMISSIBILITY, MARKER QUALITY
%  ======================================================================
function M = marker_quality (etaK, errK, theta)
etaK = etaK(:);  errK = errK(:);
M = empty_mq ();
n = numel (etaK);
if n < 2 || all (etaK == 0) || all (errK == 0), return; end
M.rho_rank = pearson_ (rank_ (etaK), rank_ (errK));
e2 = errK.^2;  tot = sum (e2);
if tot <= 0, return; end
sel = dorfler_set (etaK, theta);
M.capture     = sum (e2(sel)) / tot;
M.frac_marked = numel (sel) / n;
M.frac_ideal  = numel (dorfler_set (errK, theta)) / n;
end
function sel = dorfler_set (v, theta)
[sv, ord] = sort (v(:), 'descend');
cs = cumsum (sv.^2);
if cs(end) <= 0, sel = []; return; end
nt = find (cs >= theta*cs(end), 1, 'first');
if isempty (nt), nt = numel (sv); end
sel = ord(1:nt);
end
function r = rank_ (v)
[~, ord] = sort (v(:));  r = zeros (numel (v),1);  r(ord) = 1:numel (v);
end
function c = pearson_ (x, y)
x = x(:)-mean (x(:));  y = y(:)-mean (y(:));
d = sqrt (sum (x.^2)*sum (y.^2));
if d <= 0, c = NaN; else, c = sum (x.*y)/d; end
end

function ratio = capture_ratio (ind_K, err_K, theta)
%CAPTURE_RATIO  Error energy in the Doerfler set / energy of the optimal same-size set.
e2 = err_K(:).^2;  tot = sum (e2);
if tot <= 0, ratio = NaN; return; end
sel = dorfler_set (ind_K, theta);
nmk = max (numel (sel), 1);
cap = sum (e2(sel)) / tot;
es  = sort (e2, 'descend');
opt = sum (es(1:min (nmk, numel (es)))) / tot;
ratio = cap / max (opt, realmin);
end

function marked = mark_elements (hmsh, est, lev_of, pos_of, prm)
marked = cell (hmsh.nlevels, 1);
for l = 1:hmsh.nlevels, marked{l} = []; end
if strcmpi (prm.refinement, 'uniform')
  for l = 1:hmsh.nlevels, marked{l} = hmsh.active{l}(:).'; end
  return;
end
if all (est == 0), return; end
switch lower (prm.mark_strategy)
  case 'maximum',  sel = find (est >= prm.theta*max (est));
  case 'dorfler',  sel = dorfler_set (est, prm.theta);     % paper (4.15)
  otherwise, error ('iga:mark','Unknown mark_strategy ''%s''.', prm.mark_strategy);
end
for l = 1:hmsh.nlevels
  marked{l} = unique (pos_of(sel(lev_of(sel) == l))).';
end
end

function mf = mark_functions (hmsh, hspace, eta_K, prm)
%MARK_FUNCTIONS  eta_i^2 = sum_{K in supp(N_i)} eta_K^2, then Doerfler on functions.
%   Function marking deepens hierarchy gradually; element marking births whole levels at once.
etaf = zeros (hspace.ndof, 1);
off = 0;  ndofs = 0;
for ilev = 1:hmsh.nlevels
  ndofs = ndofs + hspace.ndof_per_level(ilev);
  nel_lev = hmsh.nel_per_level(ilev);
  if nel_lev == 0, continue; end
  sp_c = eval_conn (hmsh, hspace, ilev, hmsh.active{ilev}(:));
  C = hspace.Csub{ilev};
  conn = sp_c.connectivity;
  nel  = size (conn, 2);
  rows = conn(conn > 0);
  cols = repmat (1:nel, size (conn,1), 1);
  cols = cols(conn > 0);
  S = sparse (rows, cols, 1, size (C,1), nel);
  W = double ((C.' * S) ~= 0);
  nrw = size (W, 1);
  etaf(1:nrw) = etaf(1:nrw) + W * (eta_K(off+(1:nel)).^2);
  off = off + nel_lev;
end
mf = cell (hmsh.nlevels, 1);
for l = 1:hmsh.nlevels, mf{l} = []; end
if all (etaf == 0), return; end
[se, ord] = sort (etaf, 'descend');
cs = cumsum (se);
ntk = find (cs >= prm.theta*cs(end), 1, 'first');
if isempty (ntk), ntk = numel (se); end
sel = ord(1:ntk);
lo = 0;
for l = 1:hspace.nlevels
  nl = hspace.ndof_per_level(l);
  loc = sel(sel > lo & sel <= lo+nl) - lo;
  if ~isempty (loc)
    act = hspace.active{l}(:).';
    loc = loc(loc >= 1 & loc <= numel (act));
    mf{l} = unique (act(loc));
  end
  lo = lo + nl;
end
end

function marked = enforce_admissibility (hmsh, hspace, marked, prm)
for lev = hmsh.nlevels:-1:2
  if isempty (marked{lev}), continue; end
  ltar = lev - prm.adm_class + 1;
  if ltar < 1 || isempty (hmsh.active{ltar}), continue; end
  nd_l = nel_dir_of_level (hmsh, hspace, lev);
  nd_t = nel_dir_of_level (hmsh, hspace, ltar);
  r = round (nd_l ./ nd_t);
  if any (r < 1), continue; end
  if isempty (prm.adm_q),        qq = prm.degree + 1;
  elseif prm.adm_q == 0,         qq = max (1, ceil ((prm.degree+1)/max (r(1),1)));
  else,                          qq = prm.adm_q;
  end
  [il, jl] = ind2sub (nd_l, marked{lev}(:));
  it = ceil (il./r(1));  jt = ceil (jl./r(2));
  cand = [];
  for di = -qq:qq
    for dj = -qq:qq
      ii = min (max (it+di,1), nd_t(1));
      jj = min (max (jt+dj,1), nd_t(2));
      cand = [cand; sub2ind(nd_t, ii, jj)];   %#ok<AGROW>
    end
  end
  add = intersect (unique (cand).', hmsh.active{ltar}(:).');
  marked{ltar} = unique ([marked{ltar}(:).', add]);
end
end

function L = local_effectivity (E, etaK_in, pb, prm)
etaK = etaK_in(:);  errK = E.err_K(:);
floorv = prm.loc_tol * max (errK);
good = errK > floorv;
raw = nan (size (etaK));
raw(good) = etaK(good) ./ errK(good);
L.n_masked = sum (~good);
rc = hypot (E.cx(:) - pb.xc(1), E.cy(:) - pb.xc(2));
hd = 0.5*hypot (E.hx(:), E.hy(:));
L.is_layer = abs (rc - pb.r0) < hd + 2*pb.delta;
L.raw = raw;  L.eta = etaK;  L.err = errK;
L.cx = E.cx(:);  L.cy = E.cy(:);  L.hx = E.hx(:);  L.hy = E.hy(:);
L.med_raw_layer  = nanmedian_ (raw(L.is_layer));
L.med_raw_smooth = nanmedian_ (raw(~L.is_layer));
end
function m = nanmedian_ (v)
v = v(isfinite (v));
if isempty (v), m = NaN; else, m = median (v); end
end


%% ======================================================================
%  PROBLEM
%  ======================================================================
function pb = ex53_problem (prm)
ev = prm.eps;
pb.eps = ev;  pb.b = prm.b(:);  pb.c = prm.c;
pb.bnorm = norm (pb.b);  pb.ahat = pb.b/pb.bnorm;
pb.c0 = pb.c;
pb.r0 = prm.r0;  pb.xc = prm.xc(:);
pb.delta = sqrt (ev);
K = 2/sqrt (ev);
x0 = pb.xc(1);  y0 = pb.xc(2);  r0s = pb.r0^2;
pb.A = @(x,y) K*(r0s - (x-x0).^2 - (y-y0).^2);
pb.u = @(x,y) 16*x.*(1-x).*y.*(1-y) .* (0.5 + atan (pb.A (x,y))/pi);
pb.ug = @(x,y) ex53_grad (x, y, K, x0, y0, r0s);
pb.f = @(x,y) ex53_f (x, y, ev, pb.b, pb.c, K, x0, y0, r0s);
end
function [u, ux, uy] = ex53_grad (x, y, K, x0, y0, r0s)
a = K*(r0s - (x-x0).^2 - (y-y0).^2);  D = 1 + a.^2;
Ax = -2*K*(x-x0);  Ay = -2*K*(y-y0);
th = 0.5 + atan (a)/pi;  thx = Ax./(pi*D);  thy = Ay./(pi*D);
P = x - x.^2;  Q = y - y.^2;
ps = 16*P.*Q;  psx = 16*(1-2*x).*Q;  psy = 16*P.*(1-2*y);
u = ps.*th;  ux = psx.*th + ps.*thx;  uy = psy.*th + ps.*thy;
end
function val = ex53_f (x, y, ev, b, c, K, x0, y0, r0s)
a = K*(r0s - (x-x0).^2 - (y-y0).^2);  D = 1 + a.^2;
Ax = -2*K*(x-x0);  Ay = -2*K*(y-y0);  Axx = -2*K;  Ayy = -2*K;
th = 0.5 + atan (a)/pi;
thx = Ax./(pi*D);  thy = Ay./(pi*D);
thxx = Axx./(pi*D) - 2*a.*Ax.^2./(pi*D.^2);
thyy = Ayy./(pi*D) - 2*a.*Ay.^2./(pi*D.^2);
P = x - x.^2;  Q = y - y.^2;
ps = 16*P.*Q;  psx = 16*(1-2*x).*Q;  psy = 16*P.*(1-2*y);
psxx = -32*Q;  psyy = -32*P;
lap = psxx.*th + 2*psx.*thx + ps.*thxx + psyy.*th + 2*psy.*thy + ps.*thyy;
ux = psx.*th + ps.*thx;   uy = psy.*th + ps.*thy;
val = -ev*lap + b(1)*ux + b(2)*uy + c*(ps.*th);
end
function check_boundary_data (pb)
t = linspace (0,1,2001).';
v = [abs(pb.u(0*t,t)); abs(pb.u(0*t+1,t)); abs(pb.u(t,0*t)); abs(pb.u(t,0*t+1))];
if max (v) > 1e-12
  warning ('iga:bc','max|u_exact| on dOmega is %.3e (expected 0).', max (v));
end
end


%% ======================================================================
%  QUADRATURE HELPERS
%  ======================================================================
function tf = circle_crosses (cx, cy, hx, hy, pb, band)
rc = hypot (cx - pb.xc(1), cy - pb.xc(2));
hd = 0.5*hypot (hx, hy);
tf = abs (rc - pb.r0) < hd + band;
end
function ns = subcell_count (h, pb, prm)
ns = min (max (ceil (h/(0.5*pb.delta)), 1), prm.lay_nsub);
end
function [xq, wq] = subcell_rule (a, b, ns, nq)
[xg, wg] = gauss_ref (nq);
e = linspace (a, b, ns+1);
xq = zeros (1, ns*nq);  wq = xq;
for i = 1:ns
  z = (i-1)*nq + (1:nq);
  xq(z) = 0.5*(e(i)+e(i+1)) + 0.5*(e(i+1)-e(i))*xg;
  wq(z) = 0.5*(e(i+1)-e(i))*wg;
end
end
function [x, w] = gauss_ref (n)
persistent CA
if isempty (CA), CA = containers.Map ('KeyType','double','ValueType','any'); end
if CA.isKey (n), z = CA(n); x = z{1}; w = z{2}; return; end
k = 1:n-1;  bb = k./sqrt (4*k.^2 - 1);
T = diag (bb,1) + diag (bb,-1);
[V, D] = eig (T);
[x, ix] = sort (diag (D));
w = 2*(V(1,ix).^2);
x = x(:).';  w = w(:).';
CA(n) = {x, w};
end


%% ======================================================================
%  GEOPDES HELPERS
%  ======================================================================
function ridx = csub_row_indices (hspace, ilev, nrows)
[cri, ok] = get_prop (hspace, 'Csub_row_indices');
if ok && iscell (cri) && numel (cri) >= ilev && numel (cri{ilev}) == nrows
  ridx = cri{ilev}(:).';  return;
end
cand = [];
[act, ok] = get_prop (hspace, 'active');
if ok && iscell (act) && numel (act) >= ilev, cand = union (cand, act{ilev}(:).'); end
[dea, ok] = get_prop (hspace, 'deactivated');
if ok && iscell (dea) && numel (dea) >= ilev, cand = union (cand, dea{ilev}(:).'); end
if numel (cand) == nrows, ridx = sort (cand(:).'); return; end
ridx = [];
end
function [v, ok] = get_prop (obj, name)
v = [];  ok = false;
try, v = obj.(name);  ok = true;  catch, ok = false;  end
end
function ridx = row_map_from_connectivity (hmsh, hspace, ilev, nrows)
ridx = zeros (1, nrows);
if ~exist ('change_connectivity_localized_Csub','file'), return; end
msh_lev = msh_evaluate_element_list (mesh_of_level (hmsh, ilev), hmsh.active{ilev});
sp_raw = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', false);
conn_tp = sp_raw.connectivity;
sp_loc  = change_connectivity_localized_Csub (sp_raw, hspace, ilev);
conn_rw = sp_loc.connectivity;
g = conn_rw > 0 & conn_rw <= nrows;
ridx(conn_rw(g)) = conn_tp(g);
end
function coef = level_coefficients (hmsh, hspace, ilev, cl, ndof_tp)
nrows = numel (cl);
if nrows == ndof_tp, coef = cl(:); return; end
ridx = csub_row_indices (hspace, ilev, nrows);
if isempty (ridx), ridx = row_map_from_connectivity (hmsh, hspace, ilev, nrows); end
if isempty (ridx) || ~any (ridx > 0)
  error ('iga:csub','level %d: could not identify the tensor functions.', ilev);
end
coef = zeros (ndof_tp, 1);
m = ridx > 0;
coef(ridx(m)) = cl(m);
end
function m = mesh_of_level (hmsh, ilev)
mol = hmsh.mesh_of_level;
if iscell (mol), m = mol{ilev}; else, m = mol(ilev); end
end
function brk = level_breaks (hmsh, hspace, lev)
brk = {};
try
  m = mesh_of_level (hmsh, lev);
  if isfield (m, 'breaks'), brk = m.breaks; end
catch
  brk = {};
end
if isempty (brk)
  kn = hspace.space_of_level(lev).knots;
  brk = cell (1, numel (kn));
  for d = 1:numel (kn), brk{d} = unique (kn{d}); end
end
end
function nd = nel_dir_of_level (hmsh, hspace, lev)
brk = level_breaks (hmsh, hspace, lev);
nd = [numel(brk{1})-1, numel(brk{2})-1];
end
function [hxe, hye, cxe, cye] = elem_geom_of_level (hmsh, hspace, lev, elems)
brk = level_breaks (hmsh, hspace, lev);
bx = brk{1}(:);  by = brk{2}(:);
dx = diff (bx);  dy = diff (by);
[ix, iy] = ind2sub ([numel(dx), numel(dy)], elems(:));
hxe = reshape (dx(ix), [], 1);   hye = reshape (dy(iy), [], 1);
cxe = reshape (bx(ix)+dx(ix)/2, [], 1);
cye = reshape (by(iy)+dy(iy)/2, [], 1);
end
function d = dirichlet_dofs (hmsh, hspace)
d = [];
try
  for iside = 1:2*hmsh.ndim
    d = union (d, hspace.boundary(iside).dofs(:).');
  end
catch
  d = [];
end
if ~isempty (d), d = d(:).'; return; end
ndofs = 0;
for lev = 1:hspace.nlevels
  sp_lev = hspace.space_of_level(lev);
  bd = [];
  for iside = 1:2*hmsh.ndim
    bd = union (bd, sp_lev.boundary(iside).dofs(:).');
  end
  [~, ia] = intersect (hspace.active{lev}(:).', bd);
  d = union (d, ndofs + ia(:).');
  ndofs = ndofs + numel (hspace.active{lev});
end
d = d(:).';
end
function [N, Gx, Gy, Lp] = shape_blocks (sp_lev, iel, nqn, nsh)
N  = reshape (sp_lev.shape_functions(:,:,iel), nqn, nsh);
Gx = reshape (sp_lev.shape_function_gradients(1,:,:,iel), nqn, nsh);
Gy = reshape (sp_lev.shape_function_gradients(2,:,:,iel), nqn, nsh);
Lp = reshape (sp_lev.shape_function_laplacians(:,:,iel), nqn, nsh);
end

%% ======================================================================
%  TABLES
%  ======================================================================
function print_tables (R, prm)
%PRINT_TABLES  iter/dof/eta/||u-u_h||/I_eff per (degree, driver); 'vms'=(4.13), 'gradind'=(4.16).
ok = find ([R.ok]);
if isempty (ok), return; end
sep = repmat ('-', 1, 68);
tno = 0;
for k = ok
  tno = tno + 1;
  [nd, ev, iv] = driving_curve (R(k));
  er = R(k).err(:);  n = numel (nd);
  bigI = max (iv(isfinite (iv))) >= 1e3;
  fprintf ('\n  %s\n', sep);
  fprintf ('  Table %d.  %s,  p = %d (C^%d),  eps = %.0e,  %s refinement.\n', ...
           tno, driver_title (R(k).driver), R(k).degree, R(k).regularity, ...
           prm.eps, lower (prm.refinement));
  fprintf ('  %s\n', sep);
  fprintf ('  %4s  %8s  %13s  %15s  %11s\n', ...
           'k', 'dof', eta_head (R(k).driver), '||u-u_h||_L2', 'I_eff');
  fprintf ('  %s\n', sep);
  for it = 1:n
    if bigI
      fprintf ('  %4d  %8d  %13.4e  %15.4e  %11.3e\n', ...
               it, nd(it), ev(it), er(it), iv(it));
    else
      fprintf ('  %4d  %8d  %13.4e  %15.4e  %11.4f\n', ...
               it, nd(it), ev(it), er(it), iv(it));
    end
  end
  fprintf ('  %s\n', sep);
  rate = -fit_slope_ (sqrt (nd), er);
  if bigI
    fprintf ('  final ||u-u_h||_L2 = %.4e at %d dof;  I_eff in [%.3e, %.3e];  L2 rate vs sqrt(dof) = %.2f\n', ...
             er(end), nd(end), min (iv), max (iv), rate);
  else
    fprintf ('  final ||u-u_h||_L2 = %.4e at %d dof;  I_eff in [%.4f, %.4f];  L2 rate vs sqrt(dof) = %.2f\n', ...
             er(end), nd(end), min (iv), max (iv), rate);
  end
  nb = sum (~logical (R(k).trust(:)));
  if nb > 0
    fprintf ('  note: %d of %d iterations are at the linear-solver accuracy floor.\n', nb, n);
  end
end
fprintf ('\n');
if ~isempty (prm.tex_export), write_tex_tables (R, prm); end
end

function [nd, ev, iv] = driving_curve (Rk)
%DRIVING_CURVE  dof, driving estimator, and I_eff for this run.
nd = Rk.ndof(:);
if strcmpi (Rk.driver, 'gradind')
  ev = Rk.eta_grad(:);  iv = Rk.ieff_grad(:);
else
  ev = Rk.eta(:);       iv = Rk.ieff(:);
end
end

function s = driver_title (drv)
switch lower (drv)
  case 'gradind', s = 'gradient indicator (4.16)';
  case 'vms',     s = 'VMS residual estimator (4.13)';
  otherwise,      s = 'uniform refinement';
end
end
function s = eta_head (drv)
if strcmpi (drv,'gradind'), s = 'eta_grad'; else, s = 'eta_vms'; end
end

function write_tex_tables (R, prm)
%WRITE_TEX_TABLES  Same five columns as print_tables, written as booktabs LaTeX.
fn = prm.tex_export;
[d,~,~] = fileparts (fn);
if ~isempty (d) && ~exist (d,'dir'), mkdir (d); end
fid = fopen (fn, 'w');
if fid < 0, warning ('iga:tex','could not open %s for writing.', fn); return; end
ok = find ([R.ok]);
for k = ok
  [nd, ev, iv] = driving_curve (R(k));
  er = R(k).err(:);
  if strcmpi (R(k).driver,'gradind')
    eh = '$\eta_{\mathrm{grad}}$';  ih = '$I_{\mathrm{eff}}$';  ifmt = '%.3e';
  else
    eh = '$\eta$';                  ih = '$I_{\mathrm{eff}}$';  ifmt = '%.4f';
  end
  fprintf (fid, '\\begin{table}[htbp]\n\\centering\n');
  fprintf (fid, '\\caption{%s, $p=%d$ ($C^{%d}$), $\\varepsilon=10^{%d}$.}\n', ...
           driver_title (R(k).driver), R(k).degree, R(k).regularity, ...
           round (log10 (prm.eps)));
  fprintf (fid, '\\begin{tabular}{rrccc}\n\\toprule\n');
  fprintf (fid, '$k$ & $N$ & %s & $\\|u-u_h\\|_{L^2}$ & %s \\\\\n\\midrule\n', eh, ih);
  for it = 1:numel (nd)
    fprintf (fid, ['%d & %d & %.4e & %.4e & ' ifmt ' \\\\\n'], ...
             it, nd(it), ev(it), er(it), iv(it));
  end
  fprintf (fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');
end
fclose (fid);
fprintf ('  LaTeX tables written: %s\n', fn);
end

function s = fit_slope_ (x, y)
lx = log (x(:));  ly = log (y(:));
v = isfinite (lx) & isfinite (ly) & (x(:) > 0) & (y(:) > 0);
if sum (v) < 2, s = NaN; return; end
cf = polyfit (lx(v), ly(v), 1);  s = cf(1);
end


%% ======================================================================
%  SUMMARY STRUCT  (returned in OUT, not printed)
%  ======================================================================
function SUM = build_summary (R, prm)   %#ok<INUSD>
SUM = struct ('degree',{},'regularity',{},'driver',{},'ndof',{},'err',{}, ...
              'eta',{},'ieff',{},'ieff_min',{},'ieff_max',{},'n_below1',{}, ...
              'dof_1e4',{},'rate',{},'mq_loc',{},'mq_grad',{}, ...
              'n_excluded',{},'n_nonmono',{});
m = 0;
for k = 1:numel (R)
  if ~R(k).ok || isempty (R(k).ndof), continue; end
  m = m + 1;  e = numel (R(k).ndof);
  [~, ev, iv] = driving_curve (R(k));
  SUM(m).degree = R(k).degree;  SUM(m).regularity = R(k).regularity;
  SUM(m).driver = R(k).driver;  SUM(m).ndof = R(k).ndof(e);
  SUM(m).err = R(k).err(e);     SUM(m).eta = ev(e);
  SUM(m).ieff = iv(e);
  tr = logical (R(k).trust(:));
  if ~any (tr), tr = true (size (tr)); end
  SUM(m).ieff_min = min (iv(tr));  SUM(m).ieff_max = max (iv(tr));
  SUM(m).n_below1 = sum (iv(tr) < 1);
  SUM(m).n_excluded = sum (~logical (R(k).trust(:)));
  SUM(m).n_nonmono  = sum (diff (R(k).err(:)) > 0);
  kk = find (R(k).err < 1e-4, 1, 'first');
  if isempty (kk), SUM(m).dof_1e4 = NaN; else, SUM(m).dof_1e4 = R(k).ndof(kk); end
  SUM(m).rate = sub_first (fit_rate (R(k).ndof(:), R(k).err(:)));
  SUM(m).mq_loc = R(k).mq_loc;  SUM(m).mq_grad = R(k).mq_grad;
end
end


%% ======================================================================
%  FIGURES
%  ======================================================================
function C = degree_colors ()
C = [0.00 0.45 0.74;      % p = 1  blue
     0.85 0.33 0.10;      % p = 2  vermilion
     0.47 0.67 0.19;      % p = 3  green
     0.49 0.18 0.56;      % p = 4  purple
     0.93 0.69 0.13];     % p = 5  ochre
end
function M = degree_markers ()
M = {'o','s','^','d','v'};
end
function ls = driver_style (drv)
if strcmpi (drv,'gradind'), ls = '--'; else, ls = '-'; end
end
function s = drv_name (drv)
switch lower (drv)
  case 'vms',     s = '\eta_{vms}-driven';
  case 'gradind', s = '\eta_{grad}-driven';
  otherwise,      s = 'uniform';
end
end
function s = short_label (Rk)
switch lower (Rk.driver)
  case 'gradind', s = sprintf ('p = %d,  \\eta_{grad}', Rk.degree);
  case 'vms',     s = sprintf ('p = %d,  \\eta_{vms}', Rk.degree);
  otherwise,      s = sprintf ('p = %d,  uniform', Rk.degree);
end
end
function [degs, drvs, IDX] = run_grid (R)
ok = find ([R.ok]);
degs = unique ([R(ok).degree]);
drvs = unique_cell ({R(ok).driver});
IDX = zeros (numel (drvs), numel (degs));
for k = ok
  IDX(find (strcmpi (drvs, R(k).driver),1), find (degs == R(k).degree,1)) = k;
end
end
function c = unique_cell (c)
o = {};
for i = 1:numel (c), if ~any (strcmpi (o, c{i})), o{end+1} = c{i}; end, end  %#ok<AGROW>
c = o;
end

function siam_axes (ax, prm)
%SIAM_AXES  Consistent print-ready axes style.
set (ax, 'FontName', prm.fig_font, 'FontSize', prm.fig_fsize, ...
         'LineWidth', 0.75, 'TickDir', 'out', 'Box', 'on', ...
         'XGrid', 'on', 'YGrid', 'on', 'Layer', 'top');
try, set (ax, 'GridLineStyle', ':', 'GridAlpha', 0.30); catch, end
try, set (ax, 'MinorGridLineStyle', 'none', 'XMinorGrid', 'off', ...
              'YMinorGrid', 'off'); catch, end
end

function n = cascade_counter (resetflag)
%CASCADE_COUNTER  Step counter for figure cascading; call with true to reset.
persistent step
if isempty (step), step = 0; end
if nargin > 0 && resetflag, step = 0; end
n = step;
step = step + 1;
end

function set_fig_size (fig, w, h)
%SET_FIG_SIZE  Fixed pixel size, cascaded from the previous figure, clamped on-screen.
try
  su = get (0,'Units');  set (0,'Units','pixels');
  ss = get (0,'ScreenSize');  set (0,'Units',su);
catch
  ss = [1 1 1366 768];
end
cx  = max (1,  round ((ss(3)-w)/2));
cy  = max (45, round ((ss(4)-h)/2));
off = mod (cascade_counter (), 14) * 28;
x = max (1,  min (cx + off, max (1,  ss(3)-w-1)));
y = max (45, min (cy + off, max (45, ss(4)-h-1)));
try
  set (fig, 'Units','pixels', 'Position', [x y w h]);
catch
  set (fig, 'Units','pixels', 'Position', [100 100 w h]);
end
try, set (fig, 'PaperPositionMode','auto'); catch, end
end

function wide_log_xaxis (ax, nd, pad)
%WIDE_LOG_XAXIS  Log dof axis padded and ticked at 1-2-5 per decade.
nd = nd(isfinite (nd) & nd > 0);
if isempty (nd), return; end
lo = min (nd)/pad;  hi = max (nd)*pad;
set (ax, 'XScale', 'log');
xlim (ax, [lo hi]);
t = [];
for e = floor (log10 (lo)) : ceil (log10 (hi))
  t = [t, [1 2 5]*10^e];   %#ok<AGROW>
end
t = t(t >= lo & t <= hi);
if numel (t) >= 2
  lab = cell (1, numel (t));
  for i = 1:numel (t), lab{i} = sprintf ('%g', t(i)); end
  set (ax, 'XTick', t, 'XTickLabel', lab);
end
end

function save_figure (fig, name, prm)
if isempty (prm.fig_export), return; end
if ~exist (prm.fig_export, 'dir'), mkdir (prm.fig_export); end
try                                   % the interactive toolbar must not be
  axs = findall (fig, 'Type', 'axes');% baked into the exported figure
  for a = 1:numel (axs), set (axs(a).Toolbar, 'Visible', 'off'); end
catch
end
for i = 1:numel (prm.fig_fmt)
  f  = lower (prm.fig_fmt{i});
  fn = fullfile (prm.fig_export, [name '.' f]);
  try
    if exist ('exportgraphics', 'file')
      if strcmp (f, 'png'), exportgraphics (fig, fn, 'Resolution', 400);
      else,                 exportgraphics (fig, fn, 'ContentType', 'vector');
      end
    else
      switch f
        case 'pdf', print (fig, '-dpdf',   '-painters', fn);
        case 'eps', print (fig, '-depsc2', '-painters', fn);
        otherwise,  print (fig, '-dpng',   '-r400',     fn);
      end
    end
    fprintf ('  figure written: %s\n', fn);
  catch werr
    warning ('iga:fig', 'could not write %s (%s).', fn, werr.message);
  end
end
end

function put_legend (ax, h, leg, key, ncol)
%PUT_LEGEND  Legend below axes, column-sorted by (degree, driver).
[~, ord] = sort (key);
lg = legend (ax, h(ord), leg(ord), 'Location', 'southoutside');
try, set (lg, 'NumColumns', max (1, ncol)); catch, end
try, set (lg, 'Box', 'off', 'FontSize', 10, 'Interpreter', 'tex'); catch, end
end

function s = eps_str (ev)
e = log10 (ev);
if abs (e - round (e)) < 1e-9, s = sprintf ('10^{%d}', round (e));
else,                          s = sprintf ('%.1e', ev);
end
end

function plot_error_all (R, prm)
%PLOT_ERROR_ALL  ||u-u_h||_L2 vs dof, all degrees and drivers on one axes.
ok = find ([R.ok]);
if isempty (ok), return; end
cascade_counter (true);
fig = figure ('Name','L2 error vs dof','Color','w');
set_fig_size (fig, 980, 660);
ax = axes ('Parent', fig);  hold (ax, 'on');
C = degree_colors ();  M = degree_markers ();
h = [];  leg = {};  key = [];  ndall = [];
for k = ok
  p = min (max (R(k).degree,1), 5);
  isg = strcmpi (R(k).driver, 'gradind');
  if isg, mfc = 'w'; else, mfc = C(p,:); end
  hh = plot (ax, R(k).ndof, R(k).err, [M{p} driver_style(R(k).driver)], ...
             'Color', C(p,:), 'LineWidth', 1.6, 'MarkerSize', 6.0, ...
             'MarkerFaceColor', mfc, 'MarkerEdgeColor', C(p,:));
  h(end+1) = hh;  leg{end+1} = short_label (R(k));   %#ok<AGROW>
  key(end+1) = 10*R(k).degree + isg;                 %#ok<AGROW>
  ndall = [ndall; R(k).ndof(:)];                     %#ok<AGROW>
end
set (ax, 'XScale', 'log', 'YScale', 'log');
wide_log_xaxis (ax, ndall, prm.xpad);
xlabel (ax, 'degrees of freedom  N');
ylabel (ax, '|| u - u_h ||_{L^2(\Omega)}');
siam_axes (ax, prm);
put_legend (ax, h, leg, key, numel (unique ([R(ok).degree])));
hold (ax, 'off');
save_figure (fig, 'error_vs_dof', prm);
end

function plot_ieff_all (R, prm)
%PLOT_IEFF_ALL  I_eff of the driving estimator vs dof, all runs on one axes.
ok = find ([R.ok]);
if isempty (ok), return; end
cascade_counter (true);
fig = figure ('Name','Effectivity index vs dof','Color','w');
set_fig_size (fig, 980, 660);
ax = axes ('Parent', fig);  hold (ax, 'on');
C = degree_colors ();  M = degree_markers ();
h = [];  leg = {};  key = [];  ndall = [];
for k = ok
  p = min (max (R(k).degree,1), 5);
  [nd, ~, iv] = driving_curve (R(k));
  tr = logical (R(k).trust(:));
  if ~any (tr), tr = true (size (tr)); end
  isg = strcmpi (R(k).driver, 'gradind');
  if isg, mfc = 'w'; else, mfc = C(p,:); end
  hh = plot (ax, nd(tr), iv(tr), [M{p} driver_style(R(k).driver)], ...
             'Color', C(p,:), 'LineWidth', 1.6, 'MarkerSize', 6.0, ...
             'MarkerFaceColor', mfc, 'MarkerEdgeColor', C(p,:));
  h(end+1) = hh;  leg{end+1} = short_label (R(k));   %#ok<AGROW>
  key(end+1) = 10*R(k).degree + isg;                 %#ok<AGROW>
  ndall = [ndall; nd(tr)];                           %#ok<AGROW>
end
set (ax, 'XScale', 'log', 'YScale', 'log');
wide_log_xaxis (ax, ndall, prm.xpad);
yl  = ylim (ax);                                   % keep I_eff = 1 visible
ylo = min (0.85, 0.85*yl(1));   yhi = 1.25*yl(2);
ylim (ax, [ylo yhi]);
xl = xlim (ax);  yl = ylim (ax);
plot (ax, xl, [1 1], 'k--', 'LineWidth', 1.2);
xlim (ax, xl);  ylim (ax, yl);
xlabel (ax, 'degrees of freedom  N');
ylabel (ax, 'I_{eff} = \eta / || u - u_h ||_{L^2(\Omega)}');
siam_axes (ax, prm);
put_legend (ax, h, leg, key, numel (unique ([R(ok).degree])));
hold (ax, 'off');
save_figure (fig, 'ieff_vs_dof', prm);
end

function plot_local_maps (R, prm)
%PLOT_LOCAL_MAPS  eta_K/||u-u_h||_K on final mesh; one figure per (degree, driver).
%   Log colour scale symmetric about 1; grey = elements below the error floor.
ok = find ([R.ok]);
if isempty (ok), return; end
cm = bwr_map (255);
allv = [];
for k = ok
  v = R(k).loc.raw;  v = v(isfinite (v) & v > 0);
  allv = [allv; log10(v)];    %#ok<AGROW>
end
if isempty (allv), return; end
L = max (0.15, quantile_ (abs (allv), 0.98));
crange = [-L, L];
cascade_counter (true);
for k = ok
  fig = figure ('Name', sprintf ('Local effectivity, p=%d (%s)', ...
                                 R(k).degree, R(k).driver), 'Color','w');
  set_fig_size (fig, 560, 600);
  colormap (fig, cm);
  ax = axes ('Parent', fig);
  v = R(k).loc.raw;  v(v <= 0) = NaN;
  draw_element_map (R(k).loc, log10 (v), crange, prm);
  try, colormap (ax, cm); catch, end
  set (ax, 'FontName', prm.fig_font, 'FontSize', prm.fig_fsize);
  title (ax, sprintf ('p = %d, C^{%d},  %s mesh', R(k).degree, ...
                      R(k).regularity, drv_name (R(k).driver)));
  ratio_colorbar (ax, crange, prm);
  drawnow;
  save_figure (fig, sprintf ('local_effectivity_p%d_%s', R(k).degree, ...
                             lower (R(k).driver)), prm);
end
end

function ratio_colorbar (ax, crange, prm)
%RATIO_COLORBAR  Colour bar ticked in ratio units (1/4, 1/2, 1, 2, 4, ...).
pos = get (ax, 'Position');
cb  = colorbar ('peer', ax);
set (ax, 'Position', pos);
set (cb, 'Units', 'normalized');
set (cb, 'Position', [min(0.955, pos(1)+pos(3)+0.014), pos(2), 0.013, pos(4)]);
e  = -6:6;
tv = log10 (2.^e);
m  = tv >= crange(1) & tv <= crange(2);
tv = tv(m);  e = e(m);
lab = cell (1, numel (e));
for i = 1:numel (e)
  if e(i) >= 0, lab{i} = sprintf ('%g', 2^e(i));
  else,         lab{i} = sprintf ('1/%g', 2^(-e(i)));
  end
end
try, set (cb, 'Ticks', tv, 'TickLabels', lab); catch, end
try, set (get (cb,'Label'), 'String', '\eta_K / ||u-u_h||_K');
catch, try, ylabel (cb, '\eta_K / ||u-u_h||_K'); catch, end
end
try, set (cb, 'FontName', prm.fig_font, 'FontSize', prm.fig_fsize-2); catch, end
end

function plot_final_meshes (R, prm)
%PLOT_FINAL_MESHES  Final THB mesh with layer circle; one figure per (degree, driver).
ok = find ([R.ok]);
if isempty (ok), return; end
th = linspace (0, 2*pi, 400);
cascade_counter (true);
for k = ok
  fig = figure ('Name', sprintf ('Final mesh, p=%d (%s)', R(k).degree, ...
                                 R(k).driver), 'Color','w');
  set_fig_size (fig, 560, 600);
  ax = axes ('Parent', fig);
  [V, F] = element_boxes (R(k).loc);
  patch ('Parent',ax,'Faces',F,'Vertices',V,'FaceColor','none', ...
         'EdgeColor',[.25 .25 .25],'LineWidth',.25);
  hold (ax, 'on');
  plot (ax, prm.xc(1)+prm.r0*cos (th), prm.xc(2)+prm.r0*sin (th), ...
        'r-','LineWidth',1.4);
  axis (ax, [0 1 0 1]);  axis (ax, 'square');  box (ax, 'on');
  set (ax, 'XTick',[0 .5 1], 'YTick',[0 .5 1], 'Layer','top', ...
           'FontName', prm.fig_font, 'FontSize', prm.fig_fsize);
  title (ax, sprintf ('p = %d, C^{%d},  %s,  %d elements,  %d dof', ...
                      R(k).degree, R(k).regularity, drv_name (R(k).driver), ...
                      size (F,1), R(k).ndof(end)));
  hold (ax, 'off');
  drawnow;
  save_figure (fig, sprintf ('final_mesh_p%d_%s', R(k).degree, ...
                             lower (R(k).driver)), prm);
end
end

function [V, F, n] = element_boxes (L)
cx = L.cx(:); cy = L.cy(:); hx = L.hx(:); hy = L.hy(:);
n = numel (cx);
x0 = cx-hx/2; x1 = cx+hx/2; y0 = cy-hy/2; y1 = cy+hy/2;
V = [x0 y0; x1 y0; x1 y1; x0 y1];
F = [(1:n).', (n+1:2*n).', (2*n+1:3*n).', (3*n+1:4*n).'];
end
function draw_element_map (L, val, crange, prm)
[V, F] = element_boxes (L);
val = val(:);  bad = ~isfinite (val);
cla;
patch ('Faces',F,'Vertices',V,'FaceColor',[.86 .86 .86],'EdgeColor','none');
hold on;
if any (~bad)
  patch ('Faces',F(~bad,:),'Vertices',V,'FaceVertexCData',val(~bad), ...
         'FaceColor','flat','EdgeColor','none');
end
th = linspace (0,2*pi,400);
plot (prm.xc(1)+prm.r0*cos (th), prm.xc(2)+prm.r0*sin (th),'k-','LineWidth',1.0);
set_clim (crange);
axis ([0 1 0 1]); axis square; box on;
set (gca,'XTick',[0 .5 1],'YTick',[0 .5 1],'Layer','top');
hold off;
end
function set_clim (cr)
if ~(cr(2) > cr(1)), return; end
try, caxis (cr); catch, try, clim (cr); catch, end, end   %#ok<CAXIS>
end
function cm = bwr_map (n)
if nargin < 1, n = 255; end
m = max (2, floor (n/2));
t = linspace (0,1,m+1).';  o = ones (m+1,1);
lo = [t, t, o];  hi = [o, flipud(t), flipud(t)];
cm = [lo(1:end-1,:); hi];
end
function q = quantile_ (v, pq)
v = sort (v(isfinite (v)));
if isempty (v), q = 1; return; end
q = v(max (1, min (numel (v), round (pq*numel (v)))));
end

function plot_solution_surfs (R, prm, pb)   %#ok<INUSD>
%PLOT_SOLUTION_SURFS  u_h surface over (0,1)^2; one figure per (degree, driver).
ok = find ([R.ok]);
if isempty (ok), return; end
n = max (60, prm.plot_npts);
xg = linspace (0,1,n);
[X, Y] = ndgrid (xg, xg);
cascade_counter (true);
for k = ok
  try, Z = eval_hier_solution (R(k).sol, xg, xg); catch, continue; end
  if isempty (Z) || ~any (isfinite (Z(:))), continue; end
  fig = figure ('Name', sprintf ('Solution, p=%d (%s)', R(k).degree, ...
                                 R(k).driver), 'Color','w');
  set_fig_size (fig, 640, 600);
  ax = axes ('Parent', fig);
  surf (ax, X, Y, Z, 'EdgeColor','none');  shading (ax, 'interp');
  axis (ax, 'tight');
  xlabel (ax,'x'); ylabel (ax,'y'); zlabel (ax,'u_h');
  grid (ax,'on'); box (ax,'on'); view (ax, 135, 30);
  set (ax,'FontName',prm.fig_font,'FontSize',prm.fig_fsize);
  try, pbaspect (ax, [1 1 0.75]); catch, end
  title (ax, sprintf ('p = %d, C^{%d},  %s', R(k).degree, R(k).regularity, ...
                      drv_name (R(k).driver)));
  drawnow;
  save_figure (fig, sprintf ('solution_p%d_%s', R(k).degree, ...
                             lower (R(k).driver)), prm);
end
end

function Z = eval_hier_solution (sol, xg, yg)
hmsh = sol.hmsh;  hspace = sol.hspace;  u = sol.u;
Z = nan (numel (xg), numel (yg));
LUT = build_level_lut (hmsh, hspace);
LUT = attach_coefficients (LUT, hmsh, hspace, u);
[X, Y] = ndgrid (xg, yg);
gl = locate_hier_batch (LUT, X(:), Y(:));
lev = zeros (numel (X),1);
for l = 1:numel (LUT)
  if LUT(l).nel == 0, continue; end
  lev (gl > LUT(l).off & gl <= LUT(l).off + LUT(l).nel) = l;
end
for l = 1:numel (LUT)
  id = find (lev == l);
  if isempty (id) || isempty (LUT(l).Cm), continue; end
  V0 = tp_eval_scattered (LUT(l), X(id), Y(id));
  Z(id) = V0;
end
end