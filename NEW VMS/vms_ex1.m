function R = vms_ex1 (varargin)
%VMS_EX1  Example 1 -- 1D exponential boundary layer.  I_eff = 1.
%
%   -kappa u'' + a u' = f  on (0,1),  u(0) = u(1) = 0,
%   u = x - (e^{(x-1)/k} - e^{-1/k}) / (1 - e^{-1/k}).
%
%   ONE FILE.  No dependencies, no toolboxes, no GeoPDEs.
%
%   Edit the USER PARAMETERS block below and just run it, or override any
%   single parameter from the call without touching the file:
%   METHOD.  Solver and estimator are the SAME object -- one constrained
%   fine-scale solve on V' = { v in H^1_0 : (v', vbar') = 0 for all vbar }.
%   The solver is the exact VMS coarse problem, with no stabilisation
%   parameter,  a(ubar,vbar) + a(G'(f - L ubar), vbar) = (f,vbar),  assembled
%   column by column against ONE factorisation of the fine-scale saddle
%   system and solved directly, so ubar is exactly P_{H^1_0} u.  The estimator
%   is that same solve, eta_K = ||v||_{L2(K)}.  Since u - ubar = G'(f - L ubar)
%   is an identity, I_eff = 1.  ALGORITHM.md derives every step.

% =====================================================================
%                          USER PARAMETERS
% =====================================================================
o.p       = 1:5;        % spline degrees (C^{p-1} continuity)
o.mode    = 'both';     % 'uniform' | 'adaptive' | 'both'

% --- stopping: the loop ends at whichever of these fires FIRST -------
o.nsteps  = 20;         % max refinement iterations       (Inf to disable)
o.maxdof  = 160;        % degree-of-freedom budget        (Inf to disable)
o.tol     = 1e-11;      % stop when ||u-u_h||  < tol      (0   to disable)
o.esttol  = 0;          % stop when eta        < esttol   (0   to disable)
                        % esttol is the criterion you can use when u is
                        % unknown -- it is what the estimator is FOR
o.stall   = 0;          % stop after TWO steps that each cut the error by
                        % less than this fraction (e.g. 0.05).  0 = off
o.maxtime = 0;          % wall-clock budget per (p,mode) run, s.  0 = off

% --- mesh and adaptivity --------------------------------------------
o.nel0    = 8;          % elements in the initial mesh
o.marking = 'doerfler'; % 'doerfler' | 'maximum'
o.theta   = 0.5;        % doerfler: capture this fraction of total eta^2
                        % maximum : refine where eta_K >= theta*max(eta)
o.markmax = 1;          % never refine more than this FRACTION of the
                        % elements in one step (1 = no cap)

% --- problem data ----------------------------------------------------
o.kappa   = 1e-4;       % diffusion  (smaller = thinner layer, harder)
o.avel    = 1;          % advection velocity
o.fval    = 1;          % constant source term

% --- HOW MUCH FINE SCALE: the space V'_h the estimator/solver live in --
o.ns      = 0;          % uniform sub-cells per coarse element
                        % 0 = automatic (8; 16 at p = 1)
o.nfan    = 0;          % extra geometrically graded cells into each
                        % element's OUTFLOW node, so V'_h can hold the
                        % element sublayer.  0 = automatic (halve until the
                        % cell is kappa/|a|, at most 60 cells)
o.fsq     = 0;          % Gauss points per sub-cell.  0 = automatic (p+2)
                        % V'_h is a DISCRETISATION of an infinite-dimensional
                        % space; ns/nfan/fsq say how well it is resolved.
                        % I_eff is insensitive once the sublayer is captured,
                        % which a sweep demonstrates:
                        %   for k=4:4:32, vms_ex1('ns',k,'p',3); end

% --- reporting -------------------------------------------------------
o.rate    = 'dof';      % convergence order measured against
                        %   'dof'      err ~ dof^{-r}        (cost)
                        %   'sqrtdof'  err ~ sqrt(dof)^{-r}  (mesh size, 1/h)
o.verbose = true;       % print the per-step table
o.plot    = true;       % error / effectivity / order curves
o.plotmesh  = true;    % initial, middle and final mesh for degree PLOTP,
                        % with u behind it and the element-size distribution
o.plotspace = true;    % the two spaces: Vbar, the locally supported members
                        % of V', G'(psi_K), and the Green's function g'(x,y)
o.plotproj  = true;    % u, its H^1_0 and L^2 projections, the fine-scale
                        % field u' = u - Pu (exact vs computed), and the
                        % ELEMENTWISE certificate eta_K vs ||u-u_h||_{L2(K)}
o.plotp     = 3;        % which degree the three figures above use
o.figsave   = '';       % directory to write the figures to ('' = don't save)
o.figfmt    = 'pdf';    % 'pdf' | 'eps' (vector, for the paper) | 'png'
o.figfont   = 11;       % base font size of the figures
% =====================================================================

o = override_ (o, varargin{:});          % anything above can also be passed in
if isempty (o), if nargout, R = []; end,  return; end

den = 1 - exp (-abs (o.avel)/o.kappa);
uex = @(x) x - (exp ((x-1)*o.avel/o.kappa) - exp (-abs (o.avel)/o.kappa))/den;
modes = pickmodes_ (o.mode);

if o.verbose
  fprintf ('\n=== Example 1: 1D boundary layer, kappa = %g, a = %g ===\n', ...
           o.kappa, o.avel);
end
R = [];  k = 0;
for p = o.p(:).'
  fs = fsauto_ (o, p);                   % the fine-scale space settings
  for md = modes
    k = k + 1;
    brk = linspace (0, 1, o.nel0+1);
    dof = [];  nel = [];  err = [];  eta = [];  hist = {};
    why = 'nsteps';  t0 = tic;
    if o.verbose, header_ (o); end
    nit = min (o.nsteps, 10000);
    for it = 1:nit
      X = sv_ctx (brk, p, o.kappa, o.avel);
      if it > 1 && X.ndof > o.maxdof, why = 'maxdof'; break; end
      hist{end+1} = brk;                                          %#ok<AGROW>
      fh = o.fval*ones (size (X.xq));

      [uh, in1] = sv_vms_solve (X, fh, o.fval, fs);                % VMS solve
      er = sv_l2err (X, uh, uex);
      eK = sv_finescale (X, sv_resfun (X, uh, o.fval), in1.FSm, fs);   % estimate

      dof(end+1)=X.ndof; nel(end+1)=X.nel; err(end+1)=er; eta(end+1)=norm(eK); %#ok<AGROW>
      if o.verbose, row_ (p, md{1}, it, X.nel, X.ndof, err, eta, o.rate, dof); end
      [stop, why2] = stopnow_ (o, err, eta, toc (t0));
      if stop, why = why2; break; end
      if it == nit, break; end

      if strcmpi (md{1},'uniform'), sel = 1:X.nel;
      else, sel = mark_ (eK, o.marking, o.theta, o.markmax); end
      brk = bisect_ (brk, sel);
    end
    rk = pack_ (p, md{1}, dof, nel, err, eta, hist, why, o, toc (t0));
    if k == 1, R = rk; else, R(k) = rk; end
  end
end
if o.verbose,   fprintf ('\n');  summary_ ('Example 1', R, o); end
if o.plot,      plotres_ ('Example 1', R, o); end
if o.plotmesh,  plotmesh1_ (R, o); end
if o.plotspace, plotspace1_ (o); end
if o.plotproj,  plotproj1_ (o); end
end

% ----------------------------------------------------------------------
function fs = fsauto_ (o, p)
%FSAUTO_  Resolve the "how much fine scale" options for one degree.
fs.ns = o.ns;   if fs.ns == 0, fs.ns = 8;  if p == 1, fs.ns = 16; end, end
fs.nq = o.fsq;  if fs.nq == 0, fs.nq = p + 2; end
fs.nfan = o.nfan;                      % 0 = automatic, see SV_SUBBRK
fs.dmin = [];                          % filled in by SV_FINESCALE
end

function m = pickmodes_ (s)
switch lower (s)
  case 'uniform',  m = {'uniform'};
  case 'adaptive', m = {'adaptive'};
  case 'both',     m = {'uniform','adaptive'};
  otherwise, error ('vms:mode','mode must be uniform, adaptive or both');
end
end

% ======================================================================
%  SHARED: option override, marking, stopping, rates, reporting, plots
% ======================================================================
function o = override_ (o, varargin)
%OVERRIDE_  Let the caller change any USER PARAMETER without editing the file.
%   Accepts name/value pairs or one struct.  'help' prints the current values.
if numel (varargin) == 1 && ischar (varargin{1}) && ...
   any (strcmpi (varargin{1}, {'help','-h','?'}))
  fn = fieldnames (o);
  fprintf ('\ncurrent USER PARAMETERS (edit the block at the top of the file,\n');
  fprintf ('or pass name/value pairs, or one struct):\n\n');
  for i = 1:numel(fn)
    v = o.(fn{i});
    if ischar (v),        s = ['''' v ''''];
    elseif islogical (v), s = mat2str (v);
    elseif isscalar (v),  s = num2str (v);
    else,                 s = mat2str (v);
    end
    fprintf ('   %-10s = %s\n', fn{i}, s);
  end
  fprintf ('\nstopping   : whichever of nsteps / maxdof / tol / esttol /\n');
  fprintf ('             stall / maxtime fires first (0 or Inf disables one).\n');
  fprintf ('fine scale : ns, nfan, fsq control the space V''_h.\n');
  fprintf ('rate       : ''dof'' for err ~ dof^-r, ''sqrtdof'' for err ~ h^r.\n\n');
  o = [];  return;
end
if numel (varargin) == 1 && isstruct (varargin{1})
  s = varargin{1};  fn = fieldnames (s);
  for i = 1:numel(fn), o = setfield_ (o, fn{i}, s.(fn{i})); end
  return;
end
if mod (numel (varargin), 2)
  error ('vms:opts', 'options must come in name/value pairs (or one struct)');
end
for i = 1:2:numel(varargin)
  o = setfield_ (o, varargin{i}, varargin{i+1});
end
end

function o = setfield_ (o, name, val)
if ~isfield (o, name)
  error ('vms:opts', ['unknown parameter ''%s''.  Call with ''help'' to ' ...
                      'list them.'], name);
end
o.(name) = val;
end

function sel = mark_ (e, strategy, theta, cap)
%MARK_  Which entries to refine.
%   'doerfler'  smallest set carrying THETA of the total estimated error^2
%   'maximum'   every entry with eta >= THETA * max(eta)
%   CAP (optional, in (0,1]) caps the fraction of entries that may be marked.
if nargin < 4 || isempty (cap), cap = 1; end
e = e(:);
switch lower (strategy)
  case 'doerfler'
    [es, ix] = sort (e, 'descend');
    cs = cumsum (es.^2);
    sel = sort (ix(1:find (cs >= theta*cs(end), 1, 'first'))).';
  case 'maximum'
    sel = find (e >= theta*max (e)).';
  otherwise
    error ('vms:mark', 'marking must be ''doerfler'' or ''maximum''');
end
if isempty (sel), sel = find (e == max (e), 1).'; end
if cap > 0 && cap < 1 && numel (sel) > max (1, ceil (cap*numel(e)))
  nmax = max (1, ceil (cap*numel(e)));
  [~, ix] = sort (e(sel), 'descend');
  sel = sort (sel(ix(1:nmax)));
end
sel = sel(:).';
end

function b = bisect_ (b, idx)
b = sort ([b, (b(idx) + b(idx+1))/2]);
end

function [s, why] = stopnow_ (o, err, eta, el)
%STOPNOW_  Post-solve stopping tests, in priority order.
s = false;  why = '';
if o.tol > 0 && err(end) < o.tol
  s = true;  why = sprintf ('tol (||u-u_h|| < %g)', o.tol);       return
end
if o.esttol > 0 && eta(end) < o.esttol
  s = true;  why = sprintf ('esttol (eta < %g)', o.esttol);       return
end
if o.maxtime > 0 && el > o.maxtime
  s = true;  why = sprintf ('maxtime (%g s)', o.maxtime);         return
end
if o.stall > 0 && numel (err) >= 3
  q1 = err(end)/err(end-1);  q2 = err(end-1)/err(end-2);
  if q1 > 1 - o.stall && q2 > 1 - o.stall
    s = true;  why = sprintf ('stalled (< %g%% per step, twice)', 100*o.stall);
  end
end
end

function x = ratex_ (dof, base)
%RATEX_  The abscissa the convergence order is measured against.
if strcmpi (base, 'sqrtdof'), x = sqrt (dof(:).');
else,                         x = dof(:).';
end
end

function r = rates_ (dof, err, base)
%RATES_  Step-to-step observed order:  err ~ x^{-r},  x = dof or sqrt(dof).
x = ratex_ (dof, base);  err = err(:).';
r = nan (size (err));
for i = 2:numel(err)
  if err(i) > 0 && err(i-1) > 0 && x(i) > x(i-1)
    r(i) = -log (err(i)/err(i-1)) / log (x(i)/x(i-1));
  end
end
end

function rf = ratefit_ (dof, err, base, ntail)
%RATEFIT_  Least-squares order over the last NTAIL meshes (all if empty).
if nargin < 4 || isempty (ntail), ntail = numel (err); end
x = ratex_ (dof, base);  err = err(:).';
m = err > 0 & isfinite (err);
x = x(m);  err = err(m);
if numel (x) < 2, rf = NaN; return; end
j = max (1, numel(x)-ntail+1) : numel(x);
P  = polyfit (log (x(j)), log (err(j)), 1);
rf = -P(1);
end

function header_ (o)
fprintf ('\n%3s %-9s %5s %6s %6s %14s %14s %9s %8s\n', ...
         'p','mode','it','nel','dof','||u-u_h||','eta','I_eff', ...
         ratename_ (o.rate));
end

function s = ratename_ (base)
if strcmpi (base,'sqrtdof'), s = 'rate(h)'; else, s = 'rate(N)'; end
end

function row_ (p, md, it, nel, ndof, err, eta, base, dof)
r = rates_ (dof, err, base);
if isnan (r(end)), rs = '     -- '; else, rs = sprintf ('%8.2f', r(end)); end
fprintf ('%3d %-9s %5d %6d %6d %14.5e %14.5e %9.4f %s\n', ...
         p, md, it, nel, ndof, err(end), eta(end), eta(end)/err(end), rs);
end

function r = pack_ (p, md, dof, nel, err, eta, hist, why, o, el)
%PACK_  One run of the refinement loop, ready for the table and the plots.
r.p = p;  r.mode = md;  r.dof = dof;  r.nel = nel;
r.err = err;  r.eta = eta;  r.ieff = eta./err;
r.rate     = rates_   (dof, err, o.rate);
r.ratefit  = ratefit_ (dof, err, o.rate, []);
r.ratetail = ratefit_ (dof, err, o.rate, min (4, numel (err)));
r.ratebase = o.rate;  r.mesh = hist;  r.why = why;  r.time = el;
end

function summary_ (nm, R, o)
fprintf ('=== %s summary ===\n', nm);
fprintf (['settings: nel0 %d, %s(%.2f), markmax %g, fine scale ns %s / ' ...
          'nfan %s / fsq %s\n'], o.nel0, o.marking, o.theta, o.markmax, ...
          autostr_ (o.ns), autostr_ (o.nfan), autostr_ (o.fsq));
fprintf ('stopping: nsteps %g, maxdof %g, tol %g, esttol %g, stall %g, maxtime %g\n', ...
         o.nsteps, o.maxdof, o.tol, o.esttol, o.stall, o.maxtime);
if strcmpi (o.rate, 'sqrtdof')
  fprintf ('order measured against sqrt(dof) ~ 1/h\n\n');
else
  fprintf ('order measured against dof\n\n');
end
fprintf ('%3s %-9s %6s %7s %13s %13s %9s %8s %8s %7s  %s\n', ...
         'p','mode','steps','dof','||u-u_h||','eta','I_eff', ...
         'rate','rate-4','time', 'I_eff range / stop');
for j = 1:numel(R)
  r = R(j);
  if isempty (r.dof), continue; end
  fprintf ('%3d %-9s %6d %7d %13.4e %13.4e %9.4f %8s %8s %6.1fs  [%.4f, %.4f] %s\n', ...
           r.p, r.mode, numel (r.dof), r.dof(end), r.err(end), r.eta(end), ...
           r.ieff(end), ratestr_ (r.ratefit), ratestr_ (r.ratetail), r.time, ...
           min (r.ieff), max (r.ieff), r.why);
end
fprintf ('\nrate  = least-squares order over all meshes, rate-4 = over the last 4\n');
fprintf ('I_eff = eta / ||u-u_h||;  1 means the estimator IS the error\n\n');
end

function s = autostr_ (v)
if v == 0, s = 'auto'; else, s = num2str (v); end
end

function s = ratestr_ (r)
if isnan (r), s = '      --'; else, s = sprintf ('%8.2f', r); end
end


function s = sty_ (r)
if strcmpi (r.mode,'adaptive'), s = '--s'; else, s = '-o'; end
end

function j = pickrun_ (R, p, prefer)
%PICKRUN_  Index of the run with degree P (prefer the given mode).
j = [];
for i = 1:numel(R)
  if R(i).p == p && ~isempty (R(i).dof)
    if strcmpi (R(i).mode, prefer), j = i; return; end
    if isempty (j), j = i; end
  end
end
end

function idx = threemesh_ (n)
%THREEMESH_  Initial, middle and final index of a refinement history.
if n >= 3,     idx = [1, ceil(n/2), n];
elseif n == 2, idx = [1, 1, 2];
else,          idx = [1 1 1];
end
end

% ======================================================================
%  SUPPORT.  Everything below is the numerical toolkit, inlined so that this
%  file stands alone.  You should not need to touch any of it.
% ======================================================================

function [U, nel, brk] = sv_knots (brk, p, reg)
%SV_KNOTS  Open knot vector for degree P and interior regularity REG.
%
%   [U, NEL] = SV_KNOTS (BRK, P, REG) builds the open knot vector on the
%   breakpoints BRK (a strictly increasing row vector) for B-splines of
%   degree P whose interior smoothness is C^REG.  The multiplicity of every
%   interior knot is P - REG.  REG defaults to P-1 (maximal continuity).
%
%   REG = -1 gives the discontinuous space S^p_{-1}.

if nargin < 3 || isempty (reg), reg = p - 1; end
brk = brk(:).';
if any (diff (brk) <= 0)
  error ('sv_knots:brk', 'BRK must be strictly increasing.');
end
mult = p - reg;                      % interior knot multiplicity
if mult < 1 || mult > p + 1
  error ('sv_knots:reg', 'REG must satisfy -1 <= REG <= P-1.');
end
nel = numel (brk) - 1;
U = [repmat(brk(1), 1, p+1), ...
     reshape(repmat(brk(2:end-1), mult, 1), 1, []), ...
     repmat(brk(end), 1, p+1)];
end

function s = sv_findspan (n, p, u, U)
%SV_FINDSPAN  Knot span index (Piegl & Tiller A2.1), vectorised.
%   N = numel(U) - p - 2 is the index of the last basis function.  Handles
%   repeated interior knots: the span of a breakpoint is its LAST occurrence.
u  = u(:).';
br = unique (U);
ib = discretize (u, br);
ib(u <= br(1))   = 1;
ib(u >= br(end)) = numel (br) - 1;
sb = zeros (1, numel(br)-1);
for e = 1:numel(br)-1
  sb(e) = sum (U <= br(e)) - 1;          % 0-based index of the last knot = br(e)
end
s = sb(ib);
s = min (max (s, p), n);
end

function ders = sv_ders (i, u, p, nd, U)
%SV_DERS  Derivatives of the P+1 nonzero B-splines on span I (Piegl A2.3).
%   DERS is (ND+1) x (P+1):  DERS(k+1,j+1) = d^k/du^k N_{i-p+j,p}(u).

ndu = zeros (p+1, p+1);
left = zeros (1, p+1);  right = zeros (1, p+1);
ndu(1,1) = 1;
for j = 1:p
  left(j+1)  = u - U(i+2-j);
  right(j+1) = U(i+1+j) - u;
  saved = 0;
  for r = 0:j-1
    ndu(j+1,r+1) = right(r+2) + left(j-r+1);
    temp = ndu(r+1,j) / ndu(j+1,r+1);
    ndu(r+1,j+1) = saved + right(r+2)*temp;
    saved = left(j-r+1)*temp;
  end
  ndu(j+1,j+1) = saved;
end
ders = zeros (nd+1, p+1);
ders(1,:) = ndu(:,p+1).';
a = zeros (2, p+1);
for r = 0:p
  s1 = 0; s2 = 1;  a(1,1) = 1;
  for k = 1:nd
    d = 0;  rk = r - k;  pk = p - k;
    if r >= k
      a(s2+1,1) = a(s1+1,1) / ndu(pk+2,rk+1);
      d = a(s2+1,1) * ndu(rk+1,pk+1);
    end
    if rk >= -1, j1 = 1; else, j1 = -rk; end
    if r-1 <= pk, j2 = k-1; else, j2 = p-r; end
    for j = j1:j2
      a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j)) / ndu(pk+2,rk+j+1);
      d = d + a(s2+1,j+1) * ndu(rk+j+1,pk+1);
    end
    if r <= pk
      a(s2+1,k+1) = -a(s1+1,k) / ndu(pk+2,r+1);
      d = d + a(s2+1,k+1) * ndu(r+1,pk+1);
    end
    ders(k+1,r+1) = d;
    t = s1; s1 = s2; s2 = t;
  end
end
r = p;
for k = 1:nd
  ders(k+1,:) = ders(k+1,:) * r;
  r = r * (p - k);
end
end

function B = sv_evalmat (U, p, x, nd)
%SV_EVALMAT  Sparse matrices of B-spline values and derivatives at points X.
%   B{k+1}(m,j) = d^k/dx^k N_{j,p}(X(m)),  k = 0..ND.

if nargin < 4, nd = 0; end
x    = x(:);
ndof = numel (U) - p - 1;
np   = numel (x);
spn  = sv_findspan (ndof-1, p, x, U);
I = repmat ((1:np).', 1, p+1);
J = spn(:) - p + (0:p) + 1;
V = zeros (np, p+1, nd+1);
% group points by span: the recursion only needs the local knots
[su, ~, ig] = unique (spn(:));
for t = 1:numel(su)
  idx = find (ig == t);
  i = su(t);
  for m = idx.'
    d = sv_ders (i, x(m), p, nd, U);
    V(m,:,:) = d.';
  end
end
B = cell (1, nd+1);
for k = 0:nd
  B{k+1} = sparse (I(:), J(:), reshape (V(:,:,k+1), [], 1), np, ndof);
end
end

function [x, w] = sv_gauss (n)
%SV_GAUSS  N-point Gauss-Legendre rule on [-1,1] (Golub-Welsch).
k = (1:n-1).';
b = k ./ sqrt (4*k.^2 - 1);
J = diag (b, 1) + diag (b, -1);
[V, D] = eig (J);
[x, ix] = sort (diag (D));
w = 2 * (V(1,ix).').^2;
x = x(:).';  w = w(:).';
end

function [xq, wq, eidx] = sv_quad (brk, nq, nsub)
%SV_QUAD  Gauss-Legendre points and weights on the elements of BRK.
%
%   [XQ, WQ, EIDX] = SV_QUAD (BRK, NQ, NSUB) uses NQ points on each of NSUB
%   equal sub-cells of every element.  EIDX(m) is the element that carries
%   quadrature point XQ(m).

if nargin < 3 || isempty (nsub), nsub = 1; end
[xg, wg] = sv_gauss (nq);
nel = numel (brk) - 1;
xq = zeros (nel*nsub*nq, 1);  wq = xq;  eidx = xq;
c = 0;
for e = 1:nel
  a = brk(e);  b = brk(e+1);
  edge = linspace (a, b, nsub+1);
  for s = 1:nsub
    a1 = edge(s);  b1 = edge(s+1);
    jm = 0.5*(b1-a1);
    idx = c + (1:nq);
    xq(idx)   = 0.5*(a1+b1) + jm*xg;
    wq(idx)   = jm*wg;
    eidx(idx) = e;
    c = c + nq;
  end
end
end

function [xq, wq, eidx] = sv_quadg (brk, nq, nsub, grade)
%SV_QUADG  Gauss rule on NSUB sub-cells per element, graded to the RIGHT end.
%   Sub-cell edges follow t = 1 - (1-s)^GRADE, so the last sub-cell of each
%   element has length h*(1/NSUB)^GRADE -- enough to integrate an exponential
%   outflow layer without resolving it with the mesh.
if nargin < 4 || isempty (grade), grade = 1; end
[xg, wg] = sv_gauss (nq);
s = linspace (0, 1, nsub+1);
t = 1 - (1 - s).^grade;
nel = numel (brk) - 1;
xq = zeros (nel*nsub*nq, 1);  wq = xq;  eidx = xq;
c = 0;
for e = 1:nel
  a = brk(e);  hh = brk(e+1) - a;
  ed = a + hh*t;
  for m = 1:nsub
    a1 = ed(m);  b1 = ed(m+1);  jm = 0.5*(b1-a1);
    idx = c + (1:nq);
    xq(idx) = 0.5*(a1+b1) + jm*xg;  wq(idx) = jm*wg;  eidx(idx) = e;
    c = c + nq;
  end
end
end

function P = sv_legendre (x, ea, eb, qmax)
%SV_LEGENDRE  L2(ea,eb)-orthonormal Legendre polynomials of degree 0..QMAX.
t = (2*x(:) - (ea+eb)) / (eb-ea);        % t in [-1,1]
n = numel (t);
P = zeros (n, qmax+1);
P(:,1) = 1;
if qmax >= 1, P(:,2) = t; end
for k = 2:qmax
  P(:,k+1) = ((2*k-1)*t.*P(:,k) - (k-1)*P(:,k-1)) / k;
end
for k = 0:qmax                            % normalise in L2(ea,eb)
  P(:,k+1) = P(:,k+1) * sqrt ((2*k+1)/(eb-ea));
end
end

function S = sv_spaceU (U, p, nq, nsub)
%SV_SPACEU  Assemble a 1D B-spline space from an explicit open knot vector.
%   See SV_SPACE for the fields.  Quadrature uses NQ Gauss points on each of
%   NSUB sub-cells of every non-empty knot span.

if nargin < 3 || isempty (nq),   nq   = p + 1;  end
if nargin < 4 || isempty (nsub), nsub = 1;      end
U   = U(:).';
brk = unique (U);
[xq, wq, eidx] = sv_quad (brk, nq, nsub);
B = sv_evalmat (U, p, xq, 2);

S.U = U;  S.p = p;  S.brk = brk;  S.nel = numel (brk) - 1;
S.ndof = numel (U) - p - 1;
S.xq = xq;  S.wq = wq;  S.eidx = eidx;
S.B0 = B{1};  S.B1 = B{2};  S.B2 = B{3};
W = spdiags (wq, 0, numel(wq), numel(wq));
S.M = S.B0.' * W * S.B0;
S.K = S.B1.' * W * S.B1;
S.C = S.B0.' * W * S.B1;
S.dir  = [1, S.ndof];
S.free = setdiff (1:S.ndof, S.dir);
% interior regularity per breakpoint (for reference)
end

function S = sv_space (brk, p, reg, nq, nsub)
%SV_SPACE  1D B-spline space S^p_reg on the breakpoints BRK.
%
%   S = SV_SPACE (BRK, P, REG, NQ, NSUB) assembles the space and its mass
%   (S.M), H^1_0 stiffness (S.K) and advection (S.C, C(i,j) = (N_j', N_i))
%   matrices with NQ Gauss points on each of NSUB sub-cells per element.
%   REG defaults to P-1 (maximal continuity).  S.free are the interior dofs.

if nargin < 3 || isempty (reg),  reg  = p - 1;  end
if nargin < 4 || isempty (nq),   nq   = p + 3;  end
if nargin < 5 || isempty (nsub), nsub = 1;      end
U = sv_knots (brk, p, reg);
S = sv_spaceU (U, p, nq, nsub);
S.regularity = reg;
end

function X = sv_ctx (brk, p, kappa, avel, nq, nsub, grade)
%SV_CTX  Everything needed on one mesh: space, graded quadrature, operators.
if nargin < 5 || isempty (nq),    nq    = 2*p + 3; end
if nargin < 6 || isempty (nsub),  nsub  = 1;      end
if nargin < 7 || isempty (grade), grade = 1;      end

U = sv_knots (brk, p, p-1);
[xq, wq, eidx] = sv_quadg (brk, nq, nsub, grade);
B = sv_evalmat (U, p, xq, 2);
X.U = U;  X.p = p;  X.brk = brk(:).';  X.nel = numel(brk)-1;
X.ndof = numel(U) - p - 1;
X.xq = xq;  X.wq = wq;  X.eidx = eidx;
X.B0 = B{1};  X.B1 = B{2};  X.B2 = B{3};
W = spdiags (wq, 0, numel(wq), numel(wq));
X.M = X.B0.'*W*X.B0;  X.K = X.B1.'*W*X.B1;  X.C = X.B0.'*W*X.B1;
X.dir = [1, X.ndof];  X.free = 2:X.ndof-1;
X.h = diff (X.brk);
X.kappa = kappa;  X.a = avel;
X.alpha = abs(avel)*X.h/(2*kappa);
% element gather matrix:  E(K,m) = 1 if quadrature point m lies in K
X.E = sparse (eidx, 1:numel(xq), 1, X.nel, numel(xq));
end

% ======================================================================
%  THE FINE SCALE.  Sub-mesh, constraint, one factorisation, and the two
%  things built on it: the estimator and the exact VMS solver.
% ======================================================================

function brkf = sv_subbrk (brk, nsub, dmin, adir, nfan)
%SV_SUBBRK  Sub-mesh of BRK: NSUB uniform cells per element, then a geometric
%   fan into the downstream node (ADIR = +1 right, -1 left).
%   NFAN > 0 uses exactly NFAN fan cells; NFAN = 0 halves until the smallest
%   cell is <= DMIN ~ kappa/|a|, which resolves the element outflow sublayer
%   and keeps the fine-scale solve oscillation free.
if nargin < 4 || isempty (adir), adir = 1; end
if nargin < 5 || isempty (nfan), nfan = 0; end
brk = brk(:).';
brkf = brk(1);
for e = 1:numel(brk)-1
  a = brk(e);  b = brk(e+1);  hh = b - a;
  ed = a + hh*(1:nsub)/nsub;
  if adir > 0, node = b; else, node = a; end
  d = hh/nsub;  extra = [];
  if nfan > 0
    for j = 1:nfan, d = d/2;  extra(end+1) = node - adir*d; end   %#ok<AGROW>
  else
    while d > dmin && numel (extra) < 60
      d = d/2;  extra(end+1) = node - adir*d;                     %#ok<AGROW>
    end
  end
  seg = sort ([ed, extra]);
  seg = seg(seg > a + 1e-14 & seg <= b + 1e-14);
  brkf = [brkf, seg];                                             %#ok<AGROW>
end
brkf = unique (brkf);
end

function FSm = sv_fsctx (X, FSm, fs)
%SV_FSCTX  The discrete fine-scale space V'_h on the sub-mesh, its constraint
%   matrix B, and ONE factorisation of the saddle-point operator.
%
%   V'_h = { v in S^p_{p-1}(sub-mesh) cap H^1_0 : (v, w) = 0 for all w in W },
%   W = S^{p-2}_{p-3}(coarse mesh) = D^2 Vbar.  For p = 1, W is replaced by
%   point values at the interior knots (v = 0 there), the classical bubble.
%
%   The whole point of caching this: the estimator and the solver use the SAME
%   V'_h, so the sub-mesh, the constraint and the factorisation are built once
%   per mesh and re-used, which roughly halves the cost of a refinement step.

if nargin < 3 || isempty (fs)
  fs = struct ('ns', 8, 'nq', X.p+2, 'nfan', 0, 'dmin', []);
end
if isempty (fs.dmin), fs.dmin = X.kappa / max (abs (X.a), eps); end
if ~isempty (FSm) && isfield (FSm,'KKd') && isequal (FSm.brk, X.brk) && ...
   isequal (FSm.fs, fs) && FSm.p == X.p
  return                                    % nothing changed -- re-use it
end
p = X.p;
brkf = sv_subbrk (X.brk, fs.ns, fs.dmin, sign (X.a), fs.nfan);
if isfield (fs,'extra') && ~isempty (fs.extra)   % extra sub-mesh breakpoints
  brkf = unique ([brkf, fs.extra(:).']);         % (used by the figures)
end
Vf = sv_space (brkf, p, p-1, fs.nq);
% ---- the constraint  B v = 0,  i.e. v perp W ------------------------
if p >= 2
  Uw = sv_knots (X.brk, p-2, p-3);            % W = S^{p-2}_{p-3}, dim nel+p-2
  Bw = sv_evalmat (Uw, p-2, Vf.xq, 0);  Bw = Bw{1};
  Wd = spdiags (Vf.wq, 0, numel(Vf.wq), numel(Vf.wq));
  B  = Bw.' * Wd * Vf.B0;                     % B(m,j) = (N_j, w_m)
else
  Bc = sv_evalmat (Vf.U, p, X.brk(2:end-1).', 0);
  B  = Bc{1};                                 % v(interior knots) = 0
end
fr = Vf.free;
Bf = B(:, fr);
[~, Rq, eo] = qr (full (Bf).', 'vector');     % drop dependent rows
dg = abs (diag (Rq));
ncf = sum (dg > max (size (Bf))*eps*max (dg));
Bf = Bf(eo(1:ncf), :);
% ---- the saddle-point operator, factorised once ---------------------
A  = X.kappa*Vf.K + X.a*Vf.C;
KK = [A(fr,fr), Bf.'; Bf, sparse(ncf,ncf)];
FSm.brk = X.brk;  FSm.p = p;  FSm.fs = fs;  FSm.Vf = Vf;
FSm.B = B;  FSm.Bf = Bf;  FSm.ncf = ncf;  FSm.KK = KK;
FSm.KKd = decomposition (KK);
FSm.Bc  = sv_evalmat (X.U, p, Vf.xq, 2);      % coarse basis on the sub-mesh
FSm.E   = sparse (discretize (Vf.xq, X.brk), 1:numel(Vf.xq), 1, ...
                  X.nel, numel(Vf.xq));
end

function [v, lam] = sv_fssolve (FSm, F)
%SV_FSSOLVE  Apply G' to one or many assembled fine-scale loads (columns of F).
%   This is the ONLY place the fine-scale system is inverted, and it is one
%   back-substitution per column against the stored factorisation.
Vf = FSm.Vf;  fr = Vf.free;  k = size (F, 2);
sol = FSm.KKd \ [F(fr,:); zeros(FSm.ncf, k)];
v = zeros (Vf.ndof, k);  v(fr,:) = sol(1:numel(fr), :);
lam = sol(numel(fr)+1:end, :);
end

function [etaK, v, FSm] = sv_finescale (X, r, FSm, fs)
%SV_FINESCALE  Exact discrete fine-scale solve -- the VMS estimator with no
%   tau model.
%
%   Solves   v in V'_h,  a(v,w) = (r,w) for all w in V'_h,
%       V' = { v in H^1_0 : (v, wbar) = 0 for all wbar in S^{p-2}_{p-3} },
%   and returns eta_K = ||v||_{L2(K)} on the COARSE elements.
%
%   R is either a function handle x -> r(x) (the residual) or an already
%   assembled fine-scale load vector, which is what the plotting routines and
%   the solver pass in.

if nargin < 3, FSm = []; end
if nargin < 4, fs  = []; end
FSm = sv_fsctx (X, FSm, fs);
Vf  = FSm.Vf;
if isa (r, 'function_handle'), F = Vf.B0.' * (Vf.wq .* r (Vf.xq));
else,                          F = r;
end
v  = sv_fssolve (FSm, F);
uq = Vf.B0 * v;
etaK = sqrt (full (FSm.E * (Vf.wq .* uq.^2)));
end

function rf = sv_resfun (X, uh, fval)
%SV_RESFUN  Handle x -> f(x) - (-kappa u_h'' + a u_h')(x) for a coarse u_h.
rf = @(x) local (x, X, uh, fval);
end
function r = local (x, X, uh, fval)
B = sv_evalmat (X.U, X.p, x, 2);
r = fval*ones (numel(x),1) - (-X.kappa*(B{3}*uh) + X.a*(B{2}*uh));
end

function [e, eK] = sv_l2err (X, uh, uex, nq, nsub, grade)
%SV_L2ERR  ||u - u_h||_{L2} with a rule graded into the outflow layer.
if nargin < 4 || isempty (nq),    nq    = 10; end
if nargin < 5 || isempty (nsub),  nsub  = 40; end
if nargin < 6 || isempty (grade), grade = 5;  end
[xq, wq, eidx] = sv_quadg (X.brk, nq, nsub, grade);
Bc = sv_evalmat (X.U, X.p, xq, 0);
d  = Bc{1}*uh - uex (xq);
E  = sparse (eidx, 1:numel(xq), 1, X.nel, numel(xq));
eK = sqrt (full (E * (wq .* d.^2)));
e  = norm (eK);
end

function [uh, info] = sv_vms_solve (X, fh, fval, fs, FSm)
%SV_VMS_SOLVE  The exact VMS coarse problem -- no stabilisation parameter.
%
%   u = ubar + u',  u' = G'(f - L ubar), so testing the exact equation with a
%   coarse weight gives, with no modelling of any kind,
%
%       a(ubar, vbar) + a( G'(f - L ubar), vbar ) = (f, vbar)   for all vbar,
%
%   i.e. (A + B) ubar = F - g_f  with  B_ij = a( G'(-L phi_j), phi_i )  and
%   g_f = a( G'(f), phi_i ).  Every application of G' is ONE back-substitution
%   against the factorisation cached in FSm -- the same object the estimator
%   uses -- so solver and estimator are one thing and the converged ubar is
%   exactly P_{H^1_0} u.
%
%   B is assembled COLUMN BY COLUMN (all right-hand sides at once) and the
%   small dense system is solved directly.  Iterating instead stagnates: on
%   strongly graded meshes at p = 5 the coarse operator is so ill conditioned
%   that GMRES stalls near 1e-6, which silently corrupts u_h and makes the
%   estimator look wrong when it is not.

if nargin < 4, fs  = []; end
if nargin < 5, FSm = []; end
kap = X.kappa;  av = X.a;  fr = X.free;
A = kap*X.K + av*X.C;
F = X.B0.' * (X.wq .* fh);

FSm = sv_fsctx (X, FSm, fs);              % sub-mesh + constraint + factorisation
Vf  = FSm.Vf;

gf = feed (X, FSm, fval*ones (numel (Vf.xq), 1));   % response to f
b  = F(fr) - gf(fr);

n = numel (fr);                            % one column of B per coarse dof
E = zeros (X.ndof, n);  E(sub2ind ([X.ndof n], fr(:).', 1:n)) = 1;
Lq = -kap*(FSm.Bc{3}*E) + av*(FSm.Bc{2}*E);         % (L phi_j) on the sub-mesh
G  = feed (X, FSm, -Lq);                            % a( G'(-L phi_j), phi_i )
M  = A(fr,fr) + G(fr,:);

uh = zeros (X.ndof, 1);
uh(fr) = M \ b;
info.M = M;  info.FSm = FSm;
info.resid = norm (M*uh(fr) - b) / max (norm (b), eps);
end

% ----------------------------------------------------------------------
function g = feed (X, FSm, rq)
%FEED  a( G' r, vbar ) for one or many residuals (columns of RQ), sampled at
%   the sub-mesh quadrature points.  a(v,vbar) = kappa (v',vbar') + a (v',vbar).
Vf = FSm.Vf;
F  = Vf.B0.' * (Vf.wq .* rq);
v  = sv_fssolve (FSm, F);
vqx = Vf.B1 * v;
g = FSm.Bc{2}.' * (X.kappa * (Vf.wq .* vqx)) ...
  + FSm.Bc{1}.' * (X.a     * (Vf.wq .* vqx));
end
% ======================================================================
%  FIGURES.  Everything here is written for a paper: a common style, panel
%  tags (a), (b), ..., LaTeX labels, and optional vector export via FIGSAVE.
%
%    plotres_     convergence:  error, effectivity, observed order
%    plotmesh1_   the mesh history: initial, middle, final
%    plotspace1_  the two spaces: Vbar, V', G'(psi_K), g'(x,y)
%    plotproj1_   the projection: u, Pu, u', and the LOCAL certificate
% ======================================================================

function f = fig_ (nm, w, h, o)
%FIG_  A figure with the paper style applied to everything drawn in it.
f = figure ('Name', nm, 'Color', 'w', 'Units', 'pixels', ...
            'Position', [60 60 w h], 'PaperPositionMode', 'auto');
set (f, 'DefaultAxesFontName', 'Times New Roman', ...
        'DefaultAxesFontSize', o.figfont, ...
        'DefaultAxesLineWidth', 0.75, ...
        'DefaultAxesBox', 'on', ...
        'DefaultAxesTickDir', 'out', ...
        'DefaultAxesTickLabelInterpreter', 'latex', ...
        'DefaultTextInterpreter', 'latex', ...
        'DefaultLegendInterpreter', 'latex', ...
        'DefaultLineLineWidth', 1.3);
end

function ax = pax_ (m, n, k, tag)
%PAX_  One styled panel, with its (a)/(b)/... tag.
ax = subplot (m, n, k);
hold (ax, 'on');  grid (ax, 'on');
set (ax, 'GridAlpha', 0.10, 'Layer', 'top');
if nargin >= 4 && ~isempty (tag)
  text (ax, -0.19, 1.13, ['\bf(' tag ')'], 'Units', 'normalized', ...
        'Interpreter', 'tex', 'FontName', 'Times New Roman', ...
        'FontSize', get (ax, 'FontSize') + 1);
end
end

function C = pal_ ()
%PAL_  A print-safe, colour-blind-friendly palette.
C = [0.00 0.35 0.70;    % blue
     0.85 0.33 0.10;    % vermillion
     0.10 0.55 0.25;    % green
     0.55 0.20 0.62;    % purple
     0.90 0.62 0.00;    % amber
     0.35 0.35 0.35];   % grey
end

function savefig_ (f, o, name)
%SAVEFIG_  Write the figure if FIGSAVE names a directory.
if ~isfield (o, 'figsave') || isempty (o.figsave), return; end
if ~exist (o.figsave, 'dir'), mkdir (o.figsave); end
fn = fullfile (o.figsave, [name '.' o.figfmt]);
switch lower (o.figfmt)
  case {'pdf','eps'}, exportgraphics (f, fn, 'ContentType', 'vector');
  otherwise,          exportgraphics (f, fn, 'Resolution', 300);
end
fprintf ('   figure written to %s\n', fn);
end

function knots_ (brk, style)
%KNOTS_  Draw the element boundaries behind the data.
if nargin < 2, style = '-'; end
yl = ylim;
for i = 1:numel(brk)
  plot ([brk(i) brk(i)], yl, style, 'Color', [0.82 0.82 0.82], ...
        'LineWidth', 0.5, 'HandleVisibility', 'off');
end
ylim (yl);
uistack (findobj (gca,'Type','line','Color',[0.82 0.82 0.82]), 'bottom');
end

function x = plotgrid_ (brk, per)
%PLOTGRID_  PER points inside every cell of BRK, so that a curve resolved on
%   the sub-mesh (including the geometric fan) is drawn, not aliased.
x = [];
for e = 1:numel(brk)-1
  x = [x; linspace(brk(e), brk(e+1), per+1).'];                   %#ok<AGROW>
end
x = unique ([x; brk(:)]);
end

% ----------------------------------------------------------------------
function plotres_ (nm, R, o)
%PLOTRES_  Error, effectivity and observed order against the chosen abscissa.
f = fig_ ([nm ' convergence'], 1180, 380, o);
C = pal_ ();
if strcmpi (o.rate,'sqrtdof'), xl = '$\sqrt{\mathrm{dof}}\;\sim\;1/h$';
else,                          xl = 'degrees of freedom'; end
ie = [];
for j = 1:numel(R), ie = [ie, R(j).ieff]; end                     %#ok<AGROW>

pax_ (1,3,1,'a');
for j = 1:numel(R)
  r = R(j);  if isempty (r.dof), continue; end
  plot (ratex_ (r.dof,o.rate), r.err, sty_(r), 'Color', C(1+mod(r.p-1,6),:), ...
        'MarkerSize', 4.5, 'MarkerFaceColor','w', 'HandleVisibility','off');
end
% one legend entry per degree plus one per mode, instead of one per curve --
% ten entries would cover the data
lh = [];  lt = {};
for q = unique ([R.p])
  lh(end+1) = plot (nan, nan, '-', 'Color', C(1+mod(q-1,6),:), ...
                    'LineWidth', 1.6);                            %#ok<AGROW>
  lt{end+1} = sprintf ('$p = %d$', q);                            %#ok<AGROW>
end
if any (strcmpi ({R.mode},'uniform'))
  lh(end+1) = plot (nan, nan, '-o', 'Color', [.4 .4 .4], 'MarkerSize', 4.5, ...
                    'MarkerFaceColor','w');
  lt{end+1} = 'uniform';
end
if any (strcmpi ({R.mode},'adaptive'))
  lh(end+1) = plot (nan, nan, '--s', 'Color', [.4 .4 .4], 'MarkerSize', 4.5, ...
                    'MarkerFaceColor','w');
  lt{end+1} = 'adaptive';
end
set (gca,'XScale','log','YScale','log');
xlabel (xl);  ylabel ('$\|u-u_h\|_{L^2(\Omega)}$');
title ('error');
legend (lh, lt, 'Location','southwest', 'FontSize', o.figfont-2, ...
        'Box','off', 'NumColumns', 2);

pax_ (1,3,2,'b');
for j = 1:numel(R)
  r = R(j);  if isempty (r.dof), continue; end
  plot (ratex_ (r.dof,o.rate), r.ieff, sty_(r), 'Color', C(1+mod(r.p-1,6),:), ...
        'MarkerSize', 4.5, 'MarkerFaceColor','w');
end
yline (1, 'k:', 'LineWidth', 1.1);
set (gca,'XScale','log');
if ~isempty (ie)
  lo = min (0.995, min (ie));  hi = max (1.005, max (ie));
  ylim ([lo, hi] + 0.05*(hi-lo)*[-1 1]);
end
xlabel (xl);  ylabel ('$I_{\mathrm{eff}} = \eta\,/\,\|u-u_h\|$');
title ('effectivity');

pax_ (1,3,3,'c');
for j = 1:numel(R)
  r = R(j);  if isempty (r.dof), continue; end
  plot (ratex_ (r.dof,o.rate), r.rate, sty_(r), 'Color', C(1+mod(r.p-1,6),:), ...
        'MarkerSize', 4.5, 'MarkerFaceColor','w');
end
set (gca,'XScale','log');
xlabel (xl);
if strcmpi (o.rate,'sqrtdof'), ylabel ('observed order in $h$');
else,                          ylabel ('observed order in dof'); end
title ('convergence order');
savefig_ (f, o, 'convergence');
end

% ----------------------------------------------------------------------
function plotmesh1_ (R, o)
%PLOTMESH1_  Initial, middle and final mesh for degree PLOTP, with the exact
%   solution behind the knot lines so that one can see the mesh chase the
%   layer, and the element-size distribution on a shared logarithmic axis.
j = pickrun_ (R, o.plotp, 'adaptive');
if isempty (j)
  warning ('vms:plot','no run with p = %d to draw meshes for', o.plotp);  return
end
r = R(j);  id = threemesh_ (numel (r.mesh));
uex = uex_ (o);
f = fig_ (sprintf ('mesh history, p = %d', r.p), 1180, 620, o);
C = pal_ ();  tags = {'a','b','c','d','e','f'};
xf = unique ([linspace(0,1,2001).'; 1 - logspace(-9,-1,800).']);  % resolve the layer
for m = 1:3
  brk = r.mesh{id(m)};  nel = numel (brk)-1;
  % The layer is kappa/|a| wide, invisible against x.  Plotting against the
  % distance to the outflow boundary on a log axis shows both the layer and
  % how far into it the mesh has actually reached.
  pax_ (2,3,m,tags{m});
  s = max (1 - xf, 1e-9);
  plot (s, uex(xf), '-', 'Color', C(6,:), 'LineWidth', 1.3, 'DisplayName','$u$');
  ylim ([-0.05 1.15]);
  for i = 1:numel(brk)-1
    plot (max(1-brk(i),1e-9)*[1 1], [-0.05 1.15], '-', 'Color', C(1,:), ...
          'LineWidth', 0.55, 'HandleVisibility','off');
  end
  xline (o.kappa/abs (o.avel), '--', 'Color', C(2,:), 'LineWidth', 1.1, ...
         'HandleVisibility','off');
  set (gca, 'XScale', 'log', 'XDir', 'reverse');
  xlim ([1e-6 1]);  xlabel ('distance to the outflow, $1-x$');  ylabel ('$u$');
  title (sprintf ('%s: iteration %d, %d elements', ...
                  subsref_ ({'initial','middle','final'}, m), id(m), nel));
  ax(m) = pax_ (2,3,3+m,tags{3+m});                                %#ok<AGROW>
  stairs ([brk(1:end-1), brk(end)], [diff(brk), diff(brk(end-1:end))], ...
          '-', 'Color', C(2,:), 'LineWidth', 1.3);
  plot (0.5*(brk(1:end-1)+brk(2:end)), diff(brk), 'o', 'Color', C(2,:), ...
        'MarkerSize', 3.5, 'MarkerFaceColor', C(2,:));
  set (gca, 'YScale', 'log');
  xlim ([0 1]);  xlabel ('$x$');  ylabel ('element size $h_K$');
  title (sprintf ('$h_{\\min} = %.2e$', min (diff (brk))));
end
hmin = min (diff (r.mesh{id(3)}));  hmax = max (diff (r.mesh{id(1)}));
set (ax, 'YLim', [10^floor(log10(hmin)), 10^ceil(log10(hmax))]);
sgtitle (sprintf (['Example 1: mesh history, $p=%d$, %s, %s marking, ' ...
                   '$\\kappa = %g$'], r.p, r.mode, o.marking, o.kappa), ...
         'Interpreter','latex');
savefig_ (f, o, sprintf ('meshes_p%d', r.p));
end

% ----------------------------------------------------------------------
function plotspace1_ (o)
%PLOTSPACE1_  The two spaces the VMS decomposition splits H^1_0 into.
%
%   (a) the coarse space Vbar = S^p_{p-1} cap H^1_0 -- its B-spline basis
%   (b) LOCALLY SUPPORTED members of V' -- the fine-scale bubbles.  They
%       exist because V' is cut out by only nel+p-2 constraints, and they
%       are what a sub-mesh actually resolves.
%   (c) G'(psi_K): the fine-scale response to ONE element's residual mode.
%       For C^{p-1} splines it does not stay in that element -- this is the
%       reason an element-diagonal tau model cannot work.
%   (d) the fine-scale Green's function g'(x,y) for three source points.
p = o.plotp;
brk = linspace (0, 1, o.nel0+1);
X   = sv_ctx (brk, p, o.kappa, o.avel);
fs  = fsauto_ (o, p);
FSm = sv_fsctx (X, [], fs);
Vf  = FSm.Vf;
xf  = plotgrid_ (Vf.brk, 6);                 % resolves every sub-cell
Bf  = sv_evalmat (Vf.U, p, xf, 0);  Bf = Bf{1};
Bc  = sv_evalmat (X.U,  p, xf, 0);  Bc = Bc{1};
C   = pal_ ();
f = fig_ (sprintf ('coarse and fine scale, p = %d', p), 1150, 720, o);

% ---- (a) the coarse space -------------------------------------------
pax_ (2,2,1,'a');
for i = X.free
  plot (xf, Bc(:,i), '-', 'Color', C(1+mod(i-2,6),:), 'LineWidth', 1.1);
end
xlim ([0 1]);  knots_ (brk);
xlabel ('$x$');  ylabel ('$N_{i,p}$');
title (sprintf (['coarse space $\\bar V = S^{%d}_{%d}\\cap H^1_0$, ' ...
                 '$\\dim = %d$'], p, p-1, numel (X.free)));

% ---- (b) locally supported members of V' -----------------------------
pax_ (2,2,2,'b');
k0 = max (1, floor (X.nel/3));                       % a two-element patch
a0 = brk(k0);  b0 = brk(min(k0+2, X.nel+1));
sup = [Vf.U(1:Vf.ndof).', Vf.U((1:Vf.ndof)+p+1).'];  % support of each B-spline
in  = find (sup(:,1) >= a0-1e-12 & sup(:,2) <= b0+1e-12);
in  = intersect (in(:).', Vf.free);
Z   = null (full (FSm.B(:, in)));                    % the local part of V'
% order them by energy: the smallest Dirichlet eigenvalues of the patch,
% restricted to V'.  A raw null() basis is arbitrary and concentrates in the
% smallest sub-cells, which says nothing about the space.
Kl = Z.' * Vf.K(in,in) * Z;   Ml = Z.' * Vf.M(in,in) * Z;
[Wv, Dv] = eig (full (Kl), full (Ml), 'chol');
[~, iw]  = sort (diag (Dv), 'ascend');
nb = min (3, size (Z,2));
for m = 1:nb
  v = zeros (Vf.ndof, 1);  v(in) = Z*Wv(:,iw(m));
  y = Bf*v;  y = y / max (abs (y));
  plot (xf, y, '-', 'Color', C(m,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf ('mode %d', m));
end
xline (a0, 'k--', 'HandleVisibility','off');
xline (b0, 'k--', 'HandleVisibility','off');
xlim ([max(0,a0-0.12) min(1,b0+0.12)]);  knots_ (brk);
xlabel ('$x$');  ylabel ('normalised');
legend ('Location','southeast','FontSize',o.figfont-2,'Box','off');
title (sprintf (['fine-scale space $V''$: %d locally supported modes on ' ...
                 '2 elements'], size (Z,2)));

% ---- (c) the response to one element's residual mode -----------------
pax_ (2,2,3,'c');
ke = unique (max (1, min (X.nel, round ([0.25 0.5 0.75]*X.nel))));
Fm = zeros (Vf.ndof, numel (ke));
for m = 1:numel(ke)
  K  = ke(m);
  ix = Vf.xq >= brk(K) & Vf.xq <= brk(K+1);
  psi = zeros (size (Vf.xq));
  P = sv_legendre (Vf.xq(ix), brk(K), brk(K+1), p-1);
  psi(ix) = P(:,end);            % L2(K)-orthonormal mode of degree p-1
  Fm(:,m) = Vf.B0.' * (Vf.wq .* psi);
end
V = sv_fssolve (FSm, Fm);
for m = 1:numel(ke)
  plot (xf, Bf*V(:,m), '-', 'Color', C(m,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf ('source in $K_{%d}$', ke(m)));
end
axis tight;  yl = 1.15*max (abs (ylim))*[-1 1];  ylim (yl);
for m = 1:numel(ke)                        % shade the source element
  xr = brk(ke(m):ke(m)+1);
  patch ([xr fliplr(xr)], [yl(1) yl(1) yl(2) yl(2)], C(m,:), ...
         'FaceAlpha', 0.08, 'EdgeColor','none', 'HandleVisibility','off');
end
knots_ (brk);
xlabel ('$x$');  ylabel ('$G''(\psi_K)$');
legend ('Location','northwest','FontSize',o.figfont-2,'Box','off');
title ('$V''$ is not element-local: the response leaves its element');

% ---- (d) the fine-scale Green's function -----------------------------
pax_ (2,2,4,'d');
ys = [0.30 0.55 0.80];
gs = fs;  gs.extra = [];
dmin = X.kappa / max (abs (X.a), eps);
for m = 1:numel(ys)                       % grade the sub-mesh at the source
  d = 0.5*min (X.h);
  while d > dmin && numel (gs.extra) < 200
    d = d/2;  gs.extra = [gs.extra, ys(m)+d, ys(m)-d];
  end
end
FSg = sv_fsctx (X, [], gs);
Vg  = FSg.Vf;  xg = plotgrid_ (Vg.brk, 4);
Bg  = sv_evalmat (Vg.U, p, xg, 0);  Bg = Bg{1};
Bd  = sv_evalmat (Vg.U, p, ys(:), 0);
V   = sv_fssolve (FSg, full (Bd{1}.'));    % load of a Dirac at each y_j
for m = 1:numel(ys)
  plot (xg, Bg*V(:,m), '-', 'Color', C(m,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf ('$y = %.2f$', ys(m)));
  xline (ys(m), ':', 'Color', C(m,:), 'HandleVisibility', 'off');
end
axis tight;  knots_ (brk);
xlabel ('$x$');  ylabel ('$g''(x,y)$');
legend ('Location','northwest','FontSize',o.figfont-2,'Box','off');
title ('fine-scale Green''s function: jump at $x=y$, tail downstream');

sgtitle (sprintf (['Example 1: the two scales, $p=%d$, $\\kappa=%g$, ' ...
                   '%d elements, sub-mesh $n_s=%d$'], ...
                  p, o.kappa, X.nel, fs.ns), 'Interpreter','latex');
savefig_ (f, o, sprintf ('spaces_p%d', p));
end

% ----------------------------------------------------------------------
function plotproj1_ (o)
%PLOTPROJ1_  The projection, and the local certificate.
%
%   (a) u, its H^1_0 projection Pu, and (for contrast) its L^2 projection.
%       Pu overshoots on a coarse mesh -- that is correct, and any solver
%       that does not do it is not computing Pu.
%   (b) the fine-scale part.  The EXACT u - u_h and the COMPUTED v = G'r are
%       drawn on top of each other; the estimator reproduces the error field
%       itself, not just its norm.
%   (c) elementwise eta_K against ||u-u_h||_{L2(K)}.
%   (d) the LOCAL effectivity eta_K / ||u-u_h||_K.  Because u - ubar = G'r is
%       a pointwise identity, localising it is exact too, so this is 1 on
%       every element -- the estimator certifies the error element by element.
p = o.plotp;
brk = linspace (0, 1, o.nel0+1);
X   = sv_ctx (brk, p, o.kappa, o.avel);
fs  = fsauto_ (o, p);
[uex, duex] = uex_ (o);

fh  = o.fval*ones (size (X.xq));
[uh, in1] = sv_vms_solve (X, fh, o.fval, fs);
[etaK, v] = sv_finescale (X, sv_resfun (X, uh, o.fval), in1.FSm, fs);
[er, eK]  = sv_l2err (X, uh, uex);
ch1 = proj_ (X, uex, duex, 'h1');
cl2 = proj_ (X, uex, duex, 'l2');

Vf = in1.FSm.Vf;
xf = plotgrid_ (Vf.brk, 6);
Bf = sv_evalmat (Vf.U, p, xf, 0);  Bf = Bf{1};
Bc = sv_evalmat (X.U,  p, xf, 0);  Bc = Bc{1};
C  = pal_ ();
f = fig_ (sprintf ('projection and certificate, p = %d', p), 1150, 720, o);

% ---- (a) the projections --------------------------------------------
pax_ (2,2,1,'a');
plot (xf, uex(xf),  '-',  'Color', C(6,:), 'LineWidth', 1.8, 'DisplayName','$u$');
plot (xf, Bc*ch1,   '-',  'Color', C(1,:), 'LineWidth', 1.4, ...
      'DisplayName','$Pu$: $H^1_0$ projection');
plot (xf, Bc*cl2,   '-.', 'Color', C(3,:), 'LineWidth', 1.2, ...
      'DisplayName','$L^2$ projection');
plot (xf, Bc*uh,    ':',  'Color', C(2,:), 'LineWidth', 2.0, ...
      'DisplayName','$u_h$ (VMS solver)');
axis tight;  knots_ (brk);
xlabel ('$x$');  ylabel ('$u$');
legend ('Location','northwest','FontSize',o.figfont-2,'Box','off');
title (sprintf ('$\\|u_h - Pu\\| / \\|Pu\\| = %.1e$', ...
                norm (uh-ch1)/max (norm (ch1), eps)));

% ---- (b) the fine-scale part -----------------------------------------
pax_ (2,2,2,'b');
plot (xf, uex(xf) - Bc*uh, '-', 'Color', C(6,:), 'LineWidth', 2.2, ...
      'DisplayName','$u - u_h$ (exact)');
plot (xf, Bf*v,            '-', 'Color', C(2,:), 'LineWidth', 1.2, ...
      'DisplayName','$v = G''r$ (computed)');
axis tight;  knots_ (brk);
xlabel ('$x$');  ylabel ('$u''$');
legend ('Location','northwest','FontSize',o.figfont-2,'Box','off');
title (sprintf ('the estimator reproduces the error FIELD: $I_{\\rm eff}=%.4f$', ...
                norm (etaK)/er));

% ---- (c) elementwise ------------------------------------------------
pax_ (2,2,3,'c');
hb = bar ([eK(:), etaK(:)], 1, 'grouped');
set (hb, 'EdgeColor', 'none');
hb(1).FaceColor = C(6,:);  hb(2).FaceColor = C(2,:);
xlim ([0.4, X.nel+0.6]);
xlabel ('element $K$');  ylabel ('$L^2(K)$ norm');
legend ({'$\|u-u_h\|_{L^2(K)}$','$\eta_K$'}, 'Location','northwest', ...
        'FontSize',o.figfont-2,'Box','off');
title ('elementwise error and estimator');

% ---- (d) the local certificate ---------------------------------------
pax_ (2,2,4,'d');
loc = etaK(:)./max (eK(:), realmin);
plot (1:X.nel, loc, 'o-', 'Color', C(1,:), 'MarkerFaceColor', C(1,:), ...
      'MarkerSize', 4.5);
yline (1, 'k:', 'LineWidth', 1.1);
xlim ([0.6, X.nel+0.4]);
ylim ([min(0.99, min(loc)) max(1.01, max(loc))]);
xlabel ('element $K$');  ylabel ('$\eta_K / \|u-u_h\|_{L^2(K)}$');
title (sprintf ('local effectivity: $%.4f$ to $%.4f$', min (loc), max (loc)));

sgtitle (sprintf (['Example 1: projection and local certificate, ' ...
                   '$p=%d$, $\\kappa=%g$, %d elements'], p, o.kappa, X.nel), ...
         'Interpreter','latex');
savefig_ (f, o, sprintf ('projection_p%d', p));
end

% ----------------------------------------------------------------------
function [uex, duex] = uex_ (o)
%UEX_  The exact solution of Example 1 and its derivative.
k = o.kappa;  a = o.avel;
den = 1 - exp (-abs (a)/k);
uex  = @(x) x - (exp ((x-1)*a/k) - exp (-abs (a)/k))/den;
duex = @(x) 1 - (a/k)*exp ((x-1)*a/k)/den;
end

function c = proj_ (X, uex, duex, kind)
%PROJ_  The H^1_0 or L^2 projection of the exact solution onto Vbar.
%   A layer-graded rule is used, otherwise the projection of a 1e-4 layer is
%   itself a quadrature error.
[xq, wq] = sv_quadg (X.brk, 12, 40, 5);
B  = sv_evalmat (X.U, X.p, xq, 1);
Wd = spdiags (wq, 0, numel(wq), numel(wq));
fr = X.free;
if strcmpi (kind, 'h1')
  A = B{2}.' * (Wd * B{2});     rhs = B{2}.' * (wq .* duex (xq));
else
  A = B{1}.' * (Wd * B{1});     rhs = B{1}.' * (wq .* uex (xq));
end
c = zeros (X.ndof, 1);
c(fr) = A(fr,fr) \ rhs(fr);
end

function v = subsref_ (c, i), v = c{i}; end
