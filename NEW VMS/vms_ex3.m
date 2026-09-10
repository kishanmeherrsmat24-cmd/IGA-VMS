function R = vms_ex3 (varargin)
%VMS_EX3  Example 3 -- 2D circular interior layer (John 5.3).  I_eff = 1.
%
%   -eps Lap u + b.grad u + c u = f  on (0,1)^2,  b = (2,3), c = 2,
%   u = 16 x(1-x) y(1-y) ( 1/2 + atan(A)/pi ),
%   A = (2/sqrt eps)( 1/16 - (x-1/2)^2 - (y-1/2)^2 ),
%   an interior layer of half-width sqrt(eps) along the circle r = 1/4.
%
%   ONE FILE.  No dependencies, no toolboxes, no GeoPDEs.
%
%   Edit the USER PARAMETERS block below and just run it, or override any
%   single parameter from the call without touching the file:

%   TWO THINGS THIS EXAMPLE NEEDS THAT EXAMPLE 2 DOES NOT
%     * layer-resolving QUADRATURE (QBAND, QFAC below).  Cells whose bounding
%       box meets the annulus are subdivided until the quadrature cell is
%       about sqrt(eps)/QFAC wide.  Without it neither u nor f is integrated
%       correctly and the reported error is simply wrong -- silently.  The
%       band matters more than the cell size: at 24x24, p=3, the absolute
%       error of int u^2 is 1.0e-7 with QBAND = 5 and 6.6e-13 with QBAND = 12,
%       and 1.0e-7 is 10% of ||u-u_h||^2 once the error reaches 1e-3.
%     * a fine-scale sub-mesh able to CONTAIN the interior layer, so NS is
%       raised automatically from the mesh size and sqrt(eps).
%
%   THIS IS THE SLOWEST OF THE THREE.  Start with nel0 = 8, nsteps = 3, one
%   degree, or raise eps to 1e-3.
%
%   METHOD.  Solver and estimator are the SAME constrained fine-scale solve on
%   V' = { v in H^1_0 : (grad v, grad vbar) = 0 for all vbar }, the solver
%   being the exact VMS coarse problem with no stabilisation parameter, so
%   ubar is exactly P_{H^1_0} u and I_eff = 1.  See ALGORITHM.md.

% =====================================================================
%                          USER PARAMETERS
% =====================================================================
o.p       = 2:5;        % spline degrees (C^{p-1} continuity)
o.mode    = 'both';     % 'uniform' | 'adaptive' | 'both'

% --- stopping: the loop ends at whichever of these fires FIRST -------
o.nsteps  = 8;          % max refinement iterations       (Inf to disable)
o.maxdof  = 900;        % degree-of-freedom budget        (Inf to disable)
o.tol     = 0;          % stop when ||u-u_h||  < tol      (0   to disable)
o.esttol  = 0;          % stop when eta        < esttol   (0   to disable)
o.stall   = 0;          % stop after TWO steps that each cut the error by
                        % less than this fraction (e.g. 0.05).  0 = off
o.maxtime = 0;          % wall-clock budget per (p,mode) run, s.  0 = off
                        % worth setting here, e.g. 'maxtime', 180

% --- mesh and adaptivity --------------------------------------------
o.nel0    = 8;          % elements per direction in the initial mesh
o.marking = 'doerfler'; % 'doerfler' | 'maximum'
o.theta   = 0.35;       % doerfler: capture this fraction of total eta^2
                        % maximum : refine where the score >= theta*max
o.markdir = true;       % directional marking (score knot lines).  false marks
                        % elements and refines the union of their intervals --
                        % about 7x the dof for the same error; see VMS_EX2.
o.markmax = 1;          % never refine more than this FRACTION of the knot
                        % intervals in one step (1 = no cap).  A circular
                        % layer on a TENSOR mesh marks whole bands, so 0.3
                        % is often the better setting here.

% --- problem data ----------------------------------------------------
o.eps     = 1e-4;       % diffusion; the layer half-width is sqrt(eps)

% --- HOW MUCH FINE SCALE: the space V'_h the estimator/solver live in --
o.ns      = 0;          % uniform sub-cells per coarse element per direction
                        % 0 = automatic: enough that a sub-cell is about the
                        % layer width, ceil(h_min/sqrt(eps)), clamped to NSMAX
o.nsmax   = 12;         % cap for that automatic rule (cost grows like ns^2)
o.nfan    = 0;          % graded cells into the domain outflow.  0 = automatic
                        % (this example has no boundary layer, so it is off)
o.fsq     = 0;          % Gauss points per sub-cell.  0 = automatic (p+2)
o.substab = true;       % SUPG on the SUB-mesh; see VMS_EX2 for the numbers

% --- quadrature for the interior layer -------------------------------
o.qband   = 12;         % subdivide cells within QBAND*sqrt(eps) of r = 1/4
o.qfac    = 1;          % target quadrature cell size sqrt(eps)/QFAC
o.qmax    = 16;         % cap on the sub-cells per cell per direction

% --- reporting -------------------------------------------------------
o.rate    = 'sqrtdof';  % order against 'sqrtdof' (~1/h) or 'dof' (cost)
o.verbose = true;       % print the per-step table
o.plot    = true;       % error / effectivity / order curves
o.plotmesh  = false;    % initial, middle and final mesh for degree PLOTP
o.plotspace = false;    % the two spaces: Vbar, G'(1_K), g'(x,y_0) and u'
o.plotproj  = false;    % u, u_h = Pu, the error field u - u_h, the computed
                        % v = G'r, and the ELEMENTWISE certificate
o.plotp     = 3;        % which degree the three figures above use
o.figsave   = '';       % directory to write the figures to ('' = don't save)
o.figfmt    = 'pdf';    % 'pdf' | 'eps' (vector, for the paper) | 'png'
o.figfont   = 11;       % base font size of the figures
% =====================================================================

o = override_ (o, varargin{:});          % anything above can also be passed in
if isempty (o), if nargout, R = []; end,  return; end

D  = sv2_ex3_data (o.eps);
qs = @(ax,bx,ay,by) sv2_qsub3 (ax, bx, ay, by, D.layer, o.qband, o.qfac, o.qmax);
mk = @(bx,by,pp,nq) sv2_space_q (bx, by, pp, nq, qs);
modes = pickmodes_ (o.mode);

if o.verbose
  fprintf ('\n=== Example 3: 2D circular interior layer, eps = %g ===\n', D.eps);
  fprintf ('layer half-width sqrt(eps) = %.4g,  quadrature band %g*sqrt(eps)\n', ...
           D.layer, o.qband);
end
R = [];  k = 0;
for p = o.p(:).'
  for md = modes
    k = k + 1;
    bx = linspace (0,1,o.nel0+1);  by = bx;
    dof = [];  nel = [];  err = [];  eta = [];  hist = {};
    why = 'nsteps';  t0 = tic;
    if o.verbose, header_ (o); end
    nit = min (o.nsteps, 10000);
    FS = [];
    for it = 1:nit
      S = mk (bx, by, p, p+3);
      if it > 1 && S.ndof > o.maxdof, why = 'maxdof'; break; end
      hist{end+1} = {bx, by};                                     %#ok<AGROW>
      O = sv2_ops (S);
      fs = fsauto_ (o, p, mk, D, min (min (diff (bx)), min (diff (by))));

      [uh, in1] = sv2_vms_solve (S, O, D, FS, fs, 64);   % VMS solve
      er = sv2_l2err (S, O, uh, D);
      eK = sv2_finescale (S, uh, D, in1.FS, fs);         % VMS estimator

      dof(end+1)=S.ndof; nel(end+1)=S.nel; err(end+1)=er; eta(end+1)=norm(eK); %#ok<AGROW>
      if o.verbose, row_ (p, md{1}, it, S.nel, S.ndof, err, eta, o.rate, dof); end
      [stop, why2] = stopnow_ (o, err, eta, toc (t0));
      if stop, why = why2; break; end
      if it == nit, break; end

      if strcmpi (md{1},'uniform')
        bx = bisect_ (bx, 1:numel(bx)-1);  by = bisect_ (by, 1:numel(by)-1);
      else
        EK = reshape (eK, S.nelx, S.nely);
        if o.markdir
          bx = bisect_ (bx, mark_ (sqrt (sum (EK.^2,2)), o.marking, o.theta, o.markmax));
          by = bisect_ (by, mark_ (sqrt (sum (EK.^2,1)), o.marking, o.theta, o.markmax));
        else
          sel = mark_ (eK, o.marking, o.theta, o.markmax);
          bx = bisect_ (bx, unique (mod (sel-1, S.nelx) + 1));
          by = bisect_ (by, unique (floor ((sel-1)/S.nelx) + 1));
        end
      end
    end
    rk = pack_ (p, md{1}, dof, nel, err, eta, hist, why, o, toc (t0));
    if k == 1, R = rk; else, R(k) = rk; end
  end
end
if o.verbose,   fprintf ('\n');  summary_ ('Example 3', R, o); end
if o.plot,      plotres_ ('Example 3', R, o); end
if o.plotmesh,  plotmesh2_ ('Example 3', R, o); end
if o.plotspace
  plotspace2_ ('Example 3', D, mk, o, fsauto_ (o, o.plotp, mk, D, 1/o.nel0));
end
if o.plotproj
  plotproj2_  ('Example 3', D, mk, o, fsauto_ (o, o.plotp, mk, D, 1/o.nel0));
end
end

% ----------------------------------------------------------------------
function fs = fsauto_ (o, p, mk, D, hmin)
%FSAUTO_  Resolve the "how much fine scale" options for one degree and mesh.
%   The interior layer is sqrt(eps) wide and sits INSIDE the domain, so V'_h
%   can only contain it if a sub-cell is about that wide -- hence the
%   mesh-dependent automatic rule.
if nargin < 5 || isempty (hmin), hmin = 1/max (o.nel0,1); end
fs.ns = o.ns;
if fs.ns == 0
  fs.ns = 2;  if p == 1, fs.ns = 6; end
  fs.ns = min (o.nsmax, max (fs.ns, ceil (hmin/D.layer)));
end
fs.nq = o.fsq;  if fs.nq == 0, fs.nq = []; end       % [] = p+2, see SV2_FSCTX
fs.nfan = o.nfan;
fs.substab = o.substab;
fs.mkspace = mk;
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

function [q, w] = hb_cellquad (a, b, nq, dlay, dir, nsub)
%HB_CELLQUAD  Gauss points on [a,b], optionally NSUB sub-cells and a geometric
%   fan into the outflow boundary of the domain.
if nargin < 6 || isempty (nsub), nsub = 1; end
[g, wg] = sv_gauss (nq);
ed = a + (b-a)*(0:nsub)/nsub;
if ~isempty (dlay) && numel (dlay) >= dir && dlay(dir) > 0 && ...
   abs (b - 1) < 1e-14 && (b-a) > dlay(dir)
  d = b - a;  cuts = [];
  while d > dlay(dir) && numel (cuts) < 40
    d = d/2;  cuts(end+1) = b - d;   %#ok<AGROW>
  end
  ed = unique ([ed, cuts]);
end
q = [];  w = [];
for m = 1:numel(ed)-1
  a1 = ed(m);  b1 = ed(m+1);  jm = 0.5*(b1-a1);
  q = [q; 0.5*(a1+b1) + jm*g(:)];  w = [w; jm*wg(:)];   %#ok<AGROW>
end
end

function D = sv2_ex3_data (eps_)
%SV2_EX3_DATA  Example 3 (John 5.3): circular interior layer.
%   -eps Lap u + b.grad u + c u = f on (0,1)^2, b = (2,3), c = 2, eps = 1e-4,
%   u = 16 x(1-x) y(1-y) ( 1/2 + atan(A)/pi ),
%   A = (2/sqrt(eps)) ( 1/16 - (x-1/2)^2 - (y-1/2)^2 ),
%   an interior layer of half-width sqrt(eps) along the circle r = 1/4
%   centred at (1/2,1/2).  u vanishes on the whole boundary.
if nargin < 1 || isempty (eps_), eps_ = 1e-4; end
e = eps_;  b1 = 2;  b2 = 3;  c = 2;  s = 2/sqrt(e);
D.eps = e;  D.b = [b1; b2];  D.c = c;  D.dlay = [0 0];
D.layer = sqrt (e);
P   = @(x,y) 16*x.*(1-x).*y.*(1-y);
Px  = @(x,y) 16*(1-2*x).*y.*(1-y);
Py  = @(x,y) 16*x.*(1-x).*(1-2*y);
Pxx = @(x,y) -32*y.*(1-y);
Pyy = @(x,y) -32*x.*(1-x);
A   = @(x,y) s*(1/16 - (x-0.5).^2 - (y-0.5).^2);
Ax  = @(x,y) -2*s*(x-0.5);
Ay  = @(x,y) -2*s*(y-0.5);
Axx = -2*s;  Ayy = -2*s;
ph  = @(x,y) 0.5 + atan (A(x,y))/pi;
d   = @(x,y) 1 + A(x,y).^2;
phx = @(x,y) Ax(x,y)./(pi*d(x,y));
phy = @(x,y) Ay(x,y)./(pi*d(x,y));
phxx= @(x,y) (Axx./d(x,y) - 2*A(x,y).*Ax(x,y).^2./d(x,y).^2)/pi;
phyy= @(x,y) (Ayy./d(x,y) - 2*A(x,y).*Ay(x,y).^2./d(x,y).^2)/pi;
D.u  = @(x,y) P(x,y).*ph(x,y);
D.ux = @(x,y) Px(x,y).*ph(x,y) + P(x,y).*phx(x,y);
D.uy = @(x,y) Py(x,y).*ph(x,y) + P(x,y).*phy(x,y);
uxx  = @(x,y) Pxx(x,y).*ph(x,y) + 2*Px(x,y).*phx(x,y) + P(x,y).*phxx(x,y);
uyy  = @(x,y) Pyy(x,y).*ph(x,y) + 2*Py(x,y).*phy(x,y) + P(x,y).*phyy(x,y);
D.f  = @(x,y) -e*(uxx(x,y)+uyy(x,y)) + b1*D.ux(x,y) + b2*D.uy(x,y) + c*D.u(x,y);
D.name = 'Example 3: circular interior layer';
end

function ns = sv2_qsub3 (ax, bx, ay, by, lay, band, qfac, cap)
%SV2_QSUB3  Quadrature sub-cells for the circular interior layer of Example 3.
%   A cell is subdivided only if its bounding box meets the annulus
%   |r - 1/4| <= BAND*LAY; then the quadrature cell is brought down to about
%   LAY/QFAC, at most CAP sub-cells per direction.
%
%   BAND is the parameter that matters.  The arctan profile decays only
%   algebraically, so a band of 5 layer widths truncates a tail that is still
%   worth 1e-7 in int u^2 -- which is 10% of ||u-u_h||^2 by the time the error
%   reaches 1e-3.  Measured, 24x24 mesh, p = 3, absolute error of int u^2:
%
%       band      5          12         20
%       error   1.0e-07    6.6e-13    7.5e-13
%
%   so the default is 12.  QFAC buys nothing beyond 1 (same table, 1.4e-12 at
%   qfac = 2), because the Gauss rule already integrates arctan across a cell
%   of the layer width to machine precision.
if nargin < 6 || isempty (band), band = 12; end
if nargin < 7 || isempty (qfac), qfac = 1;  end
if nargin < 8 || isempty (cap),  cap  = 16; end
dx = max ([0.5-bx, ax-0.5, 0]);          % distance from the centre to the box
dy = max ([0.5-by, ay-0.5, 0]);
rmin = hypot (dx, dy);
rmax = max ([hypot(ax-0.5,ay-0.5), hypot(bx-0.5,ay-0.5), ...
             hypot(ax-0.5,by-0.5), hypot(bx-0.5,by-0.5)]);
lo = 0.25 - band*lay;  hi = 0.25 + band*lay;
if rmax < lo || rmin > hi
  ns = [1 1];
else
  ns = [min(cap, max (1, ceil ((bx-ax)*qfac/lay))), ...
        min(cap, max (1, ceil ((by-ay)*qfac/lay)))];
end
end

function S = sv2_space_q (brkx, brky, p, nq, qsubfun)
%SV2_SPACE_Q  Tensor space with per-cell quadrature sub-division (for the
%   interior-layer example, where the layer is not at the boundary).
% (the space is assembled field by field below -- nothing to inherit)
% rebuild the 1D rules with sub-cells wherever the layer test asks for them
for d = 1:2
  if d == 1, brk = brkx; else, brk = brky; end
  q = [];  w = [];  e = [];
  for k = 1:numel(brk)-1
    if d == 1, ns = qsubfun (brk(k), brk(k+1), 0.5, 0.5);
    else,      ns = qsubfun (0.5, 0.5, brk(k), brk(k+1)); end
    m = max (ns);
    [qq, ww] = hb_cellquad (brk(k), brk(k+1), nq, [], d, m);
    q = [q; qq];  w = [w; ww];  e = [e; k*ones(numel(qq),1)];   %#ok<AGROW>
  end
  U = sv_knots (brk, p, p-1);
  B = sv_evalmat (U, p, q, 2);
  F.U=U; F.brk=brk; F.nel=numel(brk)-1; F.ndof=numel(U)-p-1;
  F.q=q; F.w=w; F.eidx=e; F.B0=B{1}; F.B1=B{2}; F.B2=B{3};
  Wd = spdiags (w,0,numel(w),numel(w));
  F.M=F.B0.'*Wd*F.B0; F.K=F.B1.'*Wd*F.B1; F.C=F.B0.'*Wd*F.B1; F.h=diff(brk);
  if d==1, S.x=F; else, S.y=F; end
end
S.p=p; S.brkx=brkx(:).'; S.brky=brky(:).';
S.nelx=numel(brkx)-1; S.nely=numel(brky)-1; S.nel=S.nelx*S.nely;
nx=S.x.ndof; ny=S.y.ndof; S.nx=nx; S.ny=ny; S.ndof=nx*ny;
[IX,IY]=ndgrid(1:nx,1:ny); bnd=(IX==1)|(IX==nx)|(IY==1)|(IY==ny);
S.free=find(~bnd(:)).'; S.dir=find(bnd(:)).';
[EX,EY]=ndgrid(S.x.eidx,S.y.eidx); S.qelem=EX(:)+(EY(:)-1)*S.nelx;
S.Wq=kron(S.y.w,S.x.w);
S.E=sparse(S.qelem,1:numel(S.Wq),1,S.nel,numel(S.Wq));
S.hx=repmat(S.x.h(:),S.nely,1); S.hy=kron(S.y.h(:),ones(S.nelx,1));
[XQ,YQ]=ndgrid(S.x.q,S.y.q); S.xq=XQ(:); S.yq=YQ(:);
end

function O = sv2_ops (S)
%SV2_OPS  Kronecker operators and the quadrature-to-dof maps of a 2D space.
O.M = kron (S.y.M, S.x.M);
O.K = kron (S.y.M, S.x.K) + kron (S.y.K, S.x.M);
O.Cx = kron (S.y.M, S.x.C);          % (dv/dx, w) with v trial
O.Cy = kron (S.y.C, S.x.M);
% quadrature-grid basis maps (nq x ndof)
O.N   = kron (S.y.B0, S.x.B0);
O.Nx  = kron (S.y.B0, S.x.B1);
O.Ny  = kron (S.y.B1, S.x.B0);
O.Nxx = kron (S.y.B0, S.x.B2);
O.Nyy = kron (S.y.B2, S.x.B0);
end

function C = sv2_prolong1d (Uc, p, F)
%SV2_PROLONG1D  Coarse 1D B-splines expressed in a refined 1D factor F.
%   The exact operator is banded (knot insertion), so the L2 solve is
%   thresholded and returned sparse -- keeping it dense blows up the tensor
%   product in 2D.
Bcf = sv_evalmat (Uc, p, F.q, 0);  Bcf = Bcf{1};
Wq  = spdiags (F.w, 0, numel(F.w), numel(F.w));
C   = F.M \ full (F.B0.' * Wq * Bcf);
C(abs(C) < 1e-12*max(abs(C(:)))) = 0;
C = sparse (C);
res = full (max (max (abs (F.B0*C - Bcf))));
if res > 1e-9
  warning ('sv2_prolong1d:nested', 'not nested (res %.2e)', res);
end
end

function tauK = sv2_tau_supg (S, D)
%SV2_TAU_SUPG  Classical SUPG parameter on the streamline chord, h -> ell/p.
%   Works for any context that carries per-element sizes S.hx, S.hy.
b = D.b;  nb = norm (b);  bh = b/nb;  p = S.p;
ell = min (S.hx(:)/max(abs(bh(1)),1e-12), S.hy(:)/max(abs(bh(2)),1e-12));
al  = nb*ell./(2*p*D.eps);
tauK = (ell./(2*p*nb)) .* (coth (al) - 1./al);
sm = al < 1e-6;
tauK(sm) = ell(sm).*al(sm)./(6*p*nb);
end

function [e, eK] = sv2_l2err (S, O, uh, D)
%SV2_L2ERR  ||u - u_h||_{L2} on the quadrature of the context S.
d = O.N*uh - D.u (S.xq, S.yq);
eK = sqrt (full (S.E * (S.Wq .* d.^2)));
e  = norm (eK);
end

% ======================================================================
%  THE FINE SCALE IN 2D.  Sub-mesh, constraint, one factorisation, and the
%  two things built on it: the estimator and the exact VMS solver.
% ======================================================================

function FS = sv2_fsctx (S, D, FS, fs)
%SV2_FSCTX  The discrete fine-scale space V'_h on the sub-mesh, its constraint
%   matrix B, and ONE factorisation of the saddle-point operator.
%
%       V'_h = { v in V_h cap H^1_0 : (grad v, grad vbar) = 0 for all vbar },
%
%   V_h being the coarse space refined NS times per element per direction
%   (plus a geometric fan into the domain outflow, so V'_h can actually hold
%   the physical layer).  The constraint is imposed by Lagrange multipliers.
%
%   FS is cached across the estimator and the solver: they use the SAME V'_h,
%   so the sub-mesh, the prolongation, the constraint and the factorisation
%   are built once per mesh instead of twice.
%
%   FS.substab adds SUPG on the SUB-mesh.  It does not change what is being
%   computed -- solver and estimator both see it -- it only lets V'_h get away
%   with fewer cells inside the O(eps/|b|) element sublayer.  Measured on
%   Example 2, p = 3, three uniform meshes:  I_eff 1.0000 with it, 0.9998 to
%   1.0002 without, at 20% more cost.  Set fs.substab = false to check.

if isempty (fs.nq), fs.nq = S.p + 2; end
if ~isempty (FS) && isfield (FS,'KKd') && isequal (FS.brkx, S.brkx) && ...
   isequal (FS.brky, S.brky) && isequal (FS.fs, fs) && FS.p == S.p
  return                                     % nothing changed -- re-use it
end
bx = subdiv (S.brkx, fs.ns, D.dlay(1)*sign(D.b(1)), fs.nfan);
by = subdiv (S.brky, fs.ns, D.dlay(2)*sign(D.b(2)), fs.nfan);
Sf = fs.mkspace (bx, by, S.p, fs.nq);
Of = sv2_ops (Sf);
Cx = sv2_prolong1d (S.x.U, S.p, Sf.x);       % coarse basis in the fine basis
Cy = sv2_prolong1d (S.y.U, S.p, Sf.y);
Cm = kron (Cy, Cx);
B  = Cm(:, S.free).' * Of.K;                 % B(m,j) = (grad Nf_j, grad Nc_m)
FS.brkx = S.brkx;  FS.brky = S.brky;  FS.p = S.p;  FS.fs = fs;
FS.Sf = Sf;  FS.Of = Of;
FS.B  = B(:, Sf.free);
FS.Ecoarse = sparse (coarse_elem (Sf, S), 1:numel(Sf.Wq), 1, S.nel, numel(Sf.Wq));
% ---- the operator on the sub-mesh -----------------------------------
e = D.eps;  b = D.b;
A = e*Of.K + b(1)*Of.Cx + b(2)*Of.Cy + D.c*Of.M;
if fs.substab
  ts = sv2_tau_supg (Sf, D);
  tq = ts (Sf.qelem);
  P  = b(1)*Of.Nx + b(2)*Of.Ny;
  L  = -e*(Of.Nxx + Of.Nyy) + b(1)*Of.Nx + b(2)*Of.Ny + D.c*Of.N;
  Dw = spdiags (Sf.Wq .* tq, 0, numel(Sf.Wq), numel(Sf.Wq));
  A  = A + P.' * Dw * L;
  FS.Pstab = P;  FS.tq = tq;
else
  FS.Pstab = [];  FS.tq = [];
end
FS.A  = A(Sf.free, Sf.free);
FS.nc = size (FS.B, 1);
FS.KK = [FS.A, FS.B.'; FS.B, sparse(FS.nc, FS.nc)];
FS.KKd = decomposition (FS.KK);
FS.Bcf = { sv_evalmat(S.x.U, S.p, Sf.x.q, 2), sv_evalmat(S.y.U, S.p, Sf.y.q, 2) };
% coarse basis and its derivatives on the sub-mesh quadrature grid
Bx = FS.Bcf{1};  By = FS.Bcf{2};
FS.Ob.N   = kron (By{1}, Bx{1});   FS.Ob.Nx  = kron (By{1}, Bx{2});
FS.Ob.Ny  = kron (By{2}, Bx{1});   FS.Ob.Nxx = kron (By{1}, Bx{3});
FS.Ob.Nyy = kron (By{3}, Bx{1});
end

function [v, lam] = sv2_fssolve (FS, F)
%SV2_FSSOLVE  Apply G' to one or many assembled fine-scale loads (columns).
%   The ONLY place the fine-scale system is inverted: one back-substitution
%   per column against the stored factorisation.
fr = FS.Sf.free;  k = size (F, 2);
sol = FS.KKd \ [F(fr,:); zeros(FS.nc, k)];
v = zeros (FS.Sf.ndof, k);  v(fr,:) = sol(1:numel(fr), :);
lam = sol(numel(fr)+1:end, :);
end

function F = sv2_fsload (FS, rq)
%SV2_FSLOAD  (r, w) for w in the sub-mesh basis, plus the SUPG term if on.
F = FS.Of.N.' * (FS.Sf.Wq .* rq);
if ~isempty (FS.Pstab)
  F = F + FS.Pstab.' * ((FS.Sf.Wq .* FS.tq) .* rq);
end
end

function [etaK, v, FS] = sv2_finescale (S, uh, D, FS, fs)
%SV2_FINESCALE  Exact discrete fine-scale solve in 2D -- the VMS estimator
%   with no tau model.  eta_K = ||v||_{L2(K)} on the COARSE elements.
FS = sv2_fsctx (S, D, FS, fs);
Sf = FS.Sf;
rq = sv2_resq (S, FS, uh, D);
v  = sv2_fssolve (FS, sv2_fsload (FS, rq));
Vq = FS.Of.N * v;
etaK = sqrt (full (FS.Ecoarse * (Sf.Wq .* Vq.^2)));
end

function rq = sv2_resq (S, FS, uh, D)
%SV2_RESQ  The strong residual f - L u_h of a coarse solution, sampled at the
%   quadrature points of the sub-mesh.
Bx = FS.Bcf{1};  By = FS.Bcf{2};
U  = reshape (uh, S.nx, S.ny);
val = @(a,c) reshape (Bx{a}*U*By{c}.', [], 1);
Sf = FS.Sf;
uv  = val (1,1);  ux = val (2,1);  uy = val (1,2);
uxx = val (3,1);  uyy = val (1,3);
rq = D.f (Sf.xq, Sf.yq) - ( -D.eps*(uxx + uyy) + D.b(1)*ux + D.b(2)*uy + D.c*uv );
end

function brkf = subdiv (brk, nsub, dlay, nfan)
%SUBDIV  NSUB uniform sub-cells per element, plus a geometric fan into the
%   OUTFLOW boundary of the DOMAIN (not of every element -- grading every
%   element starves the interior).  NFAN > 0 forces that many fan cells.
if nargin < 4 || isempty (nfan), nfan = 0; end
brkf = brk(1);
for e = 1:numel(brk)-1
  brkf = [brkf, brk(e) + (brk(e+1)-brk(e))*(1:nsub)/nsub];        %#ok<AGROW>
end
if abs (dlay) > 0 || nfan > 0
  if dlay >= 0, x0 = brk(end); sg = 1; else, x0 = brk(1); sg = -1; end
  d = (brk(end)-brk(1))/(numel(brk)-1)/nsub;  fan = [];
  if nfan > 0
    for j = 1:nfan, d = d/2;  fan(end+1) = x0 - sg*d; end          %#ok<AGROW>
  else
    while d > abs(dlay) && numel (fan) < 45
      d = d/2;  fan(end+1) = x0 - sg*d;                            %#ok<AGROW>
    end
  end
  brkf = [brkf, fan];
end
brkf = unique (brkf);
end

function ce = coarse_elem (Sf, S)
ex = discretize (Sf.x.q, S.brkx);  ex(isnan(ex)) = S.nelx;
ey = discretize (Sf.y.q, S.brky);  ey(isnan(ey)) = S.nely;
[EX, EY] = ndgrid (ex, ey);
ce = EX(:) + (EY(:)-1)*S.nelx;
end

function [uh, info] = sv2_vms_solve (S, O, D, FS, fs, blk)
%SV2_VMS_SOLVE  The exact VMS coarse problem in 2D -- no stabilisation
%   parameter; solver and estimator become the same object.
%
%       a(ubar, vbar) + a( G'(f - L ubar), vbar ) = (f, vbar)   for all vbar,
%
%   i.e. (A + B) ubar = F - g_f.  Every application of G' is ONE constrained
%   fine-scale solve, so B is assembled column by column (a whole block of
%   right-hand sides at once against the single stored factorisation) and the
%   coarse system is solved directly.  Iterating instead stagnates on graded
%   meshes and silently corrupts u_h.
%
%   BLK is the column block size (memory / speed trade-off).

if nargin < 6 || isempty (blk), blk = 64; end
e = D.eps;  b = D.b;  c = D.c;  fr = S.free;  n = numel (fr);
A = e*O.K + b(1)*O.Cx + b(2)*O.Cy + c*O.M;
F = O.N.' * (S.Wq .* D.f (S.xq, S.yq));

FS = sv2_fsctx (S, D, FS, fs);
Sf = FS.Sf;  Ob = FS.Ob;

gf = feed2 (FS, D, D.f (Sf.xq, Sf.yq));            % response to f
rhs = F(fr) - gf(fr);

M = full (A(fr,fr));
for i0 = 1:blk:n                                    % assemble B in blocks
  idx = i0:min (n, i0+blk-1);
  E = sparse (fr(idx), 1:numel(idx), 1, S.ndof, numel(idx));
  Lq = -e*(Ob.Nxx*E + Ob.Nyy*E) + b(1)*(Ob.Nx*E) + b(2)*(Ob.Ny*E) + c*(Ob.N*E);
  G = feed2 (FS, D, -full (Lq));
  M(:, idx) = M(:, idx) + G(fr, :);
end
uh = zeros (S.ndof,1);
uh(fr) = M \ rhs;
info.FS = FS;  info.M = M;
info.resid = norm (M*uh(fr) - rhs) / max (norm (rhs), eps);
end

% ----------------------------------------------------------------------
function g = feed2 (FS, D, rq)
%FEED2  a( G' r, vbar ) for one or many residuals (columns of RQ).
%   a(v,vbar) = eps (grad v, grad vbar) + (b.grad v, vbar) + c (v,vbar).
Sf = FS.Sf;  Of = FS.Of;  Ob = FS.Ob;  Wq = Sf.Wq;
v  = sv2_fssolve (FS, sv2_fsload (FS, rq));
vx = Of.Nx*v;  vy = Of.Ny*v;  v0 = Of.N*v;
g = D.eps*(Ob.Nx.'*(Wq.*vx) + Ob.Ny.'*(Wq.*vy)) ...
    + Ob.N.'*(Wq .* (D.b(1)*vx + D.b(2)*vy)) + D.c*(Ob.N.'*(Wq.*v0));
end

% ======================================================================
%  FIGURES.  Everything here is written for a paper: a common style, panel
%  tags (a), (b), ..., LaTeX labels, and optional vector export via FIGSAVE.
%
%    plotres_     convergence:  error, effectivity, observed order
%    plotmesh2_   the mesh history: initial, middle, final
%    plotspace2_  the two spaces: Vbar, G'(1_K), g'(x,y), u'
%    plotproj2_   u, u_h, the error field, and the ELEMENTWISE certificate
% ======================================================================

function f = fig_ (nm, w, h, o)
%FIG_  A figure with the paper style applied to everything drawn in it.
f = figure ('Name', nm, 'Color', 'w', 'Units', 'pixels', ...
            'Position', [50 50 w h], 'PaperPositionMode', 'auto');
set (f, 'DefaultAxesFontName', 'Times New Roman', ...
        'DefaultAxesFontSize', o.figfont, ...
        'DefaultAxesLineWidth', 0.75, ...
        'DefaultAxesBox', 'on', ...
        'DefaultAxesTickDir', 'out', ...
        'DefaultAxesTickLabelInterpreter', 'latex', ...
        'DefaultTextInterpreter', 'latex', ...
        'DefaultLegendInterpreter', 'latex', ...
        'DefaultColorbarTickLabelInterpreter', 'latex', ...
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
C = [0.00 0.35 0.70;  0.85 0.33 0.10;  0.10 0.55 0.25;
     0.55 0.20 0.62;  0.90 0.62 0.00;  0.35 0.35 0.35];
end

function M = cdiv_ ()
%CDIV_  A blue-white-red diverging map, for SIGNED fields.
t  = linspace (0, 1, 128).';
lo = [0.02 0.19 0.42];  hi = [0.70 0.09 0.11];  w = [1 1 1];
M  = [lo + t.*(w - lo); w + t.*(hi - w)];
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

function field_ (gx, gy, Z, bx, by, ttl, sym)
%FIELD_  One field over the unit square, with the coarse mesh on top.
%   SYM = true uses a symmetric diverging map (for signed fields).
imagesc (gx, gy, Z.');  set (gca, 'YDir', 'normal');  axis square; hold on
if nargin >= 7 && sym
  m = max (abs (Z(:)));  if m == 0, m = 1; end
  clim ([-m m]);  colormap (gca, cdiv_ ());
else
  colormap (gca, parula);
end
cb = colorbar;  cb.LineWidth = 0.5;
for i = 1:numel(bx), plot ([bx(i) bx(i)],[0 1],'-','Color',[1 1 1 0.30]); end
for i = 1:numel(by), plot ([0 1],[by(i) by(i)],'-','Color',[1 1 1 0.30]); end
xlim ([0 1]);  ylim ([0 1]);  xlabel ('$x$');  ylabel ('$y$');
if nargin >= 6 && ~isempty (ttl), title (ttl); end
end

function cellmap_ (bx, by, Z, ttl, ctr)
%CELLMAP_  One value per ELEMENT, drawn on the real (non-uniform) mesh.
%   CTR, if given, centres a diverging map on that value.
Zp = zeros (numel(by), numel(bx));
Zp(1:end-1, 1:end-1) = Z.';
pcolor (bx, by, Zp);  shading flat;  axis square;  hold on
if nargin >= 5 && ~isempty (ctr)
  d = max (max (abs (Z(:) - ctr)), eps);
  clim (ctr + [-d d]);  colormap (gca, cdiv_ ());
else
  colormap (gca, parula);
end
cb = colorbar;  cb.LineWidth = 0.5;
set (gca, 'Layer', 'top');
xlim ([0 1]);  ylim ([0 1]);  xlabel ('$x$');  ylabel ('$y$');
if nargin >= 4 && ~isempty (ttl), title (ttl); end
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
xlabel (xl);  ylabel ('$\|u-u_h\|_{L^2(\Omega)}$');  title ('error');
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
function plotmesh2_ (nm, R, o)
%PLOTMESH2_  Initial, middle and final mesh for degree PLOTP, with the
%   element-size distribution of each direction underneath.
j = pickrun_ (R, o.plotp, 'adaptive');
if isempty (j)
  warning ('vms:plot','no run with p = %d to draw meshes for', o.plotp);  return
end
r = R(j);  id = threemesh_ (numel (r.mesh));
f = fig_ (sprintf ('%s meshes, p = %d', nm, r.p), 1150, 680, o);
C = pal_ ();  tags = {'a','b','c','d','e','f'};
lab = {'initial','middle','final'};
for m = 1:3
  M = r.mesh{id(m)};  bx = M{1};  by = M{2};
  pax_ (2,3,m,tags{m});  axis square
  for i = 1:numel(bx)
    plot ([bx(i) bx(i)], [0 1], '-', 'Color', C(1,:), 'LineWidth', 0.4);
  end
  for i = 1:numel(by)
    plot ([0 1], [by(i) by(i)], '-', 'Color', C(1,:), 'LineWidth', 0.4);
  end
  xlim ([0 1]); ylim ([0 1]);  xlabel ('$x$'); ylabel ('$y$');
  title (sprintf ('%s: iteration %d, $%d\\times%d$', ...
                  lab{m}, id(m), numel(bx)-1, numel(by)-1));
  ax(m) = pax_ (2,3,3+m,tags{3+m});                                %#ok<AGROW>
  plot (0.5*(bx(1:end-1)+bx(2:end)), diff(bx), 'o-', 'Color', C(2,:), ...
        'MarkerSize', 3.5, 'MarkerFaceColor', C(2,:), 'DisplayName','$h_x$');
  plot (0.5*(by(1:end-1)+by(2:end)), diff(by), 's--', 'Color', C(3,:), ...
        'MarkerSize', 3.5, 'MarkerFaceColor', C(3,:), 'DisplayName','$h_y$');
  set (gca, 'YScale', 'log');
  xlim ([0 1]);  xlabel ('$x$ or $y$');  ylabel ('element size');
  legend ('Location','southwest','FontSize',o.figfont-3,'Box','off');
  title (sprintf ('$h_{\\min} = %.2e$', min ([diff(bx), diff(by)])));
end
M0 = r.mesh{id(1)};  M1 = r.mesh{id(3)};
lo = min ([diff(M1{1}), diff(M1{2})]);  hi = max ([diff(M0{1}), diff(M0{2})]);
set (ax, 'YLim', [10^floor(log10(lo)), 10^ceil(log10(hi))]);
sgtitle (sprintf ('%s: mesh history, $p=%d$, %s, %s marking', ...
                  nm, r.p, r.mode, o.marking), 'Interpreter','latex');
savefig_ (f, o, sprintf ('meshes_p%d', r.p));
end

% ----------------------------------------------------------------------
function plotspace2_ (nm, D, mk, o, fs)
%PLOTSPACE2_  The two spaces the VMS decomposition splits H^1_0 into.
%
%   (a) one basis function of the coarse space Vbar
%   (b) G'(1_K): the fine-scale response to ONE element's residual.  The
%       advective wake is what makes V' non-local for C^{p-1} splines.
%   (c) the same for a source on a single sub-cell -- an approximation of the
%       fine-scale Green's function g'(x,y_0).
%   (d) the fine-scale field u' that the estimator integrates.
p  = o.plotp;
bx = linspace (0,1,o.nel0+1);  by = bx;
S  = mk (bx, by, p, p+3);
O  = sv2_ops (S);
FS = sv2_fsctx (S, D, [], fs);
Sf = FS.Sf;
ng = 241;  gx = linspace (0,1,ng).';  gy = gx;
Ef = { sv_evalmat(Sf.x.U, p, gx, 0), sv_evalmat(Sf.y.U, p, gy, 0) };
Ec = { sv_evalmat(S.x.U,  p, gx, 0), sv_evalmat(S.y.U,  p, gy, 0) };
ev = @(E,V,nx,ny) full (E{1}{1} * reshape (V, nx, ny) * E{2}{1}.');
f  = fig_ (sprintf ('%s: coarse and fine scale, p = %d', nm, p), 1080, 780, o);

pax_ (2,2,1,'a');
A = zeros (S.ndof,1);
A(round (S.nx/2) + (round (S.ny/2)-1)*S.nx) = 1;
field_ (gx, gy, ev (Ec, A, S.nx, S.ny), bx, by, ...
        sprintf ('coarse space: one basis function of $S^{%d}_{%d}\\otimes S^{%d}_{%d}$', ...
                 p, p-1, p, p-1), false);

pax_ (2,2,2,'b');
K = [ceil(S.nelx/3), ceil(S.nely/3)];
in = Sf.xq >= bx(K(1)) & Sf.xq <= bx(K(1)+1) & ...
     Sf.yq >= by(K(2)) & Sf.yq <= by(K(2)+1);
r1 = double (in) / max (sqrt (sum (Sf.Wq(in))), eps);
v1 = sv2_fssolve (FS, sv2_fsload (FS, r1));
field_ (gx, gy, ev (Ef, v1, Sf.nx, Sf.ny), bx, by, ...
        sprintf ('$V''$: $G''(1_K)$, $K=(%d,%d)$ -- note the wake', K(1), K(2)), true);
rectangle ('Position', [bx(K(1)) by(K(2)) diff(bx(K(1):K(1)+1)) ...
                        diff(by(K(2):K(2)+1))], 'EdgeColor','k','LineWidth',1.1);

pax_ (2,2,3,'c');
c0 = [0.35 0.35];
[~, ix] = min (abs (Sf.x.brk - c0(1)));  [~, iy] = min (abs (Sf.y.brk - c0(2)));
ax0 = Sf.x.brk(max(1,ix-1));  ax1 = Sf.x.brk(min(end,ix+1));
ay0 = Sf.y.brk(max(1,iy-1));  ay1 = Sf.y.brk(min(end,iy+1));
in = Sf.xq >= ax0 & Sf.xq <= ax1 & Sf.yq >= ay0 & Sf.yq <= ay1;
r2 = double (in) / max (sqrt (sum (Sf.Wq(in))), eps);
v2 = sv2_fssolve (FS, sv2_fsload (FS, r2));
field_ (gx, gy, ev (Ef, v2, Sf.nx, Sf.ny), bx, by, ...
        sprintf ('$g''(x,y_0)$: unit source on one sub-cell at $(%.2f,%.2f)$', ...
                 c0(1), c0(2)), true);

pax_ (2,2,4,'d');
uh = sv2_vms_solve (S, O, D, FS, fs, 64);
[etaK, v] = sv2_finescale (S, uh, D, FS, fs);
er = sv2_l2err (S, O, uh, D);
field_ (gx, gy, ev (Ef, v, Sf.nx, Sf.ny), bx, by, ...
        sprintf ('$u''$ on $%d\\times%d$: $\\eta = %.3e$, $I_{\\rm eff} = %.4f$', ...
                 S.nelx, S.nely, norm (etaK), norm (etaK)/er), true);
sgtitle (sprintf ('%s: the two scales, $p=%d$, $\\varepsilon=%g$', nm, p, D.eps), ...
         'Interpreter','latex');
savefig_ (f, o, sprintf ('spaces_p%d', p));
end

% ----------------------------------------------------------------------
function plotproj2_ (nm, D, mk, o, fs)
%PLOTPROJ2_  The projection, and the local certificate, in 2D.
%
%   (a) the exact solution
%   (b) the coarse solution u_h = Pu produced by the VMS solver
%   (c) the exact error field u - u_h
%   (d) the COMPUTED fine-scale field v = G'r -- the same picture as (c)
%   (e) the elementwise estimator eta_K
%   (f) the local effectivity eta_K / ||u-u_h||_{L2(K)}.  Because
%       u - ubar = G'r is a pointwise identity, localising it is exact, so
%       this is 1 on every element.
p  = o.plotp;
bx = linspace (0,1,o.nel0+1);  by = bx;
S  = mk (bx, by, p, p+3);
O  = sv2_ops (S);
[uh, in1] = sv2_vms_solve (S, O, D, [], fs, 64);
[etaK, v, FS] = sv2_finescale (S, uh, D, in1.FS, fs);
[er, eK] = sv2_l2err (S, O, uh, D);
Sf = FS.Sf;
ng = 241;  gx = linspace (0,1,ng).';  gy = gx;
Ef = { sv_evalmat(Sf.x.U, p, gx, 0), sv_evalmat(Sf.y.U, p, gy, 0) };
Ec = { sv_evalmat(S.x.U,  p, gx, 0), sv_evalmat(S.y.U,  p, gy, 0) };
ev = @(E,V,nx,ny) full (E{1}{1} * reshape (V, nx, ny) * E{2}{1}.');
[GX, GY] = ndgrid (gx, gy);
Uex = D.u (GX, GY);
Uh  = ev (Ec, uh, S.nx, S.ny);
f = fig_ (sprintf ('%s: projection and certificate, p = %d', nm, p), 1250, 700, o);

pax_ (2,3,1,'a');
field_ (gx, gy, Uex, bx, by, 'exact solution $u$', false);
pax_ (2,3,2,'b');
field_ (gx, gy, Uh, bx, by, ...
        sprintf ('$u_h = Pu$ on $%d\\times%d$, $p=%d$', S.nelx, S.nely, p), false);
pax_ (2,3,3,'c');
field_ (gx, gy, Uex - Uh, bx, by, ...
        sprintf ('exact error $u-u_h$, $\\|\\cdot\\| = %.3e$', er), true);
pax_ (2,3,4,'d');
field_ (gx, gy, ev (Ef, v, Sf.nx, Sf.ny), bx, by, ...
        sprintf ('computed $v = G''r$, $\\eta = %.3e$', norm (etaK)), true);
pax_ (2,3,5,'e');
cellmap_ (bx, by, reshape (etaK, S.nelx, S.nely), ...
          'elementwise estimator $\eta_K$');
pax_ (2,3,6,'f');
loc = reshape (etaK(:)./max (eK(:), realmin), S.nelx, S.nely);
cellmap_ (bx, by, loc, ...
          sprintf ('local effectivity: $%.4f$ to $%.4f$', min (loc(:)), max (loc(:))), 1);
sgtitle (sprintf (['%s: projection and local certificate, $p=%d$, ' ...
                   '$\\varepsilon=%g$, $I_{\\rm eff}=%.4f$'], ...
                  nm, p, D.eps, norm (etaK)/er), 'Interpreter','latex');
savefig_ (f, o, sprintf ('projection_p%d', p));
end

