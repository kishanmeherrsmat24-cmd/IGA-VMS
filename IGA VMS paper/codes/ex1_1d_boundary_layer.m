function out = iga_supg_vms_1d_sweep (varargin)
%IGA_SUPG_VMS_1D_SWEEP  Degree sweep p=1..5 (C^{p-1}) for
%   -kappa u'' + a u' + s u = f in (0,1),  u(0)=u(1)=0.
%   Runs UNIFORM then ADAPTIVE (Doerfler on VMS eta_K) refinement per degree.
%   Requires GeoPDEs 3.x and the NURBS toolbox.

%%                         USER PARAMETERS

prm.kappa   = 1e-4;         % diffusion
prm.a       = 1;            % advection velocity
prm.s       = 0;            % reaction
prm.f       = 1;            % source
prm.degrees = 1:5;          % which of the five fixed (p, C^{p-1}) pairs
prm.nsub0   = 14;           % initial elements
prm.nsteps_uniform  = 5;    % uniform refinement passes
prm.nsteps_adaptive = 30;   % adaptive refinement passes (upper bound)
prm.max_dof    = 130;       % adaptive stop: dof budget
prm.tol        = 1e-20;     % adaptive stop: eta tolerance
prm.test_op    = 'sgs';     % 'supg' | 'sgs' | 'gls'
prm.stabilise  = true;      % false -> plain Galerkin
prm.mark_strategy = 'dorfler';   % 'dorfler' | 'maximum'
prm.mark_param = 0.40;      % Doerfler theta
prm.nquad_stab_add = 1;     % p + 1 Gauss points for the stabilisation term
prm.nquad_add      = 5;     % p + 5 Gauss points for everything else
prm.twin_quad  = true;      % two rules on the same knots; false -> use the
                            % p+5 rule for the stabilisation term too

prm.quad_mode  = 'subcell'; % 'subcell' (layer-resolving) | 'plain'
prm.lay_bandw  = 6;         % layer band width, in units of kappa*log(1/kappa)
prm.lay_nsub   = 24;        % cap on sub-cells per element
prm.lay_grade  = 3;         % geometric grading exponent toward x = 1
prm.verbose = true;
prm.plots.error = true;
prm.plots.ieff  = true;

prm = parse_options (prm, varargin);

SWEEP = [1 0; 2 1; 3 2; 4 3; 5 4];      % [degree, regularity] - fixed pairs

prm.degrees = prm.degrees(:).';
if any (prm.degrees < 1 | prm.degrees > 5)
  error ('iga:degrees', 'prm.degrees must be a subset of 1:5.');
end
if ~any (strcmpi (prm.test_op, {'supg','sgs','gls'}))
  error ('iga:testop', 'prm.test_op must be ''supg'', ''sgs'' or ''gls''.');
end

kap = prm.kappa;
prm.uex   = @(x) x - (exp((x-1)/kap) - exp(-1/kap)) ./ (1 - exp(-1/kap));
prm.delta = kap * log (1/kap);          % boundary layer width at x = 1
if ~(prm.delta > 0), prm.delta = kap; end

modes = {'uniform', 'adaptive'};

if prm.verbose
  fprintf ('\n%s\n', repmat ('=', 1, 84));
  fprintf ([' 1D stabilised IGA degree sweep   p = 1 C^0 ... p = 5 C^4\n' ...
            ' VMS fine-scale estimator (residual-based)\n' ...
            ' -%g u'''' + %g u'' + %g u = %g on (0,1),  u(0) = u(1) = 0\n' ...
            ' every degree: (1) uniform refinement, then (2) adaptive\n' ...
            ' refinement driven by the VMS indicator -- both run automatically\n' ...
            ' test operator: ''%s''   |   theta = %.2f   |   %s quadrature\n' ...
            ' tau_stab = (h/2|a|)(coth Pe - 1/Pe)                       (all p)\n' ...
            ' tau_est  = min( h/(sqrt(2)|a|) , h^2/(3 sqrt(10) kappa) ) (all p)\n' ...
            ' eta_K = tau_est,K ||a u_h'' - kappa u_h'''' + s u_h - f||_K\n'], ...
           prm.kappa, prm.a, prm.s, prm.f, lower (prm.test_op), ...
           prm.mark_param, lower (prm.quad_mode));
  fprintf ('%s\n', repmat ('=', 1, 84));
end

%%                        THE RUNS  (degree x mode)

R = struct ('degree', {}, 'regularity', {}, 'mode', {}, 'dof', {}, ...
            'L2', {}, 'eta', {}, 'Ieff', {}, 'secs', {}, 'ok', {});

nrun = 0;
for k = 1:numel (prm.degrees)
  id  = prm.degrees(k);
  p   = SWEEP(id,1);
  reg = SWEEP(id,2);

  for j = 1:numel (modes)
    md   = modes{j};
    nrun = nrun + 1;
    R(nrun).degree = p;  R(nrun).regularity = reg;
    R(nrun).mode   = md;  R(nrun).ok = false;  R(nrun).secs = NaN;

    if strcmpi (md, 'uniform'), nsteps = prm.nsteps_uniform;
    else,                       nsteps = prm.nsteps_adaptive;
    end

    t0 = tic;
    try
      S = run_1d (p, reg, nsteps, md, prm);
      R(nrun).dof  = S.dof;
      R(nrun).L2   = S.L2;
      R(nrun).eta  = S.eta;
      R(nrun).Ieff = S.Ieff;
      R(nrun).ok   = true;
    catch err
      fprintf (2, '  p = %d, %s FAILED: %s\n', p, md, err.message);
    end
    R(nrun).secs = toc (t0);
    if prm.verbose && R(nrun).ok
      fprintf ('  [p = %d, %-9s] %6.2f s   %d dof, ||e|| = %.3e\n', ...
               p, md, R(nrun).secs, R(nrun).dof(end), R(nrun).L2(end));
    end
  end
end

%%                    PER-RUN TABLES  (iteration, dof, eta, L2, Ieff)

if prm.verbose
  for k = 1:numel (R)
    if ~R(k).ok, continue; end
    fprintf ('\n  p = %d,  C^%d,  %s refinement,  test_op = ''%s''   (%.2f s)\n', ...
             R(k).degree, R(k).regularity, R(k).mode, lower (prm.test_op), ...
             R(k).secs);
    fprintf ('  %s\n', repmat ('-', 1, 56));
    fprintf ('   it        dof      eta       L2         Ieff\n');
    fprintf ('  %s\n', repmat ('-', 1, 56));
    for it = 1:numel (R(k).dof)
      fprintf ('  %3d   %8d  %10.4e  %10.4e  %8.3f\n', ...
               it, R(k).dof(it), R(k).eta(it), R(k).L2(it), R(k).Ieff(it));
    end
    fprintf ('  %s\n', repmat ('-', 1, 56));
  end
end

print_summary (R, prm);

if nargout > 0
  out = struct ('runs', {R}, 'prm', prm, 'modes', {modes});
end

if prm.plots.error
  plot_metric (R, prm, 'uniform',  'L2', '|| u - u_h ||_{L^2(0,1)}', ...
               'L2 error vs dof (uniform)', false);
  plot_metric (R, prm, 'adaptive', 'L2', '|| u - u_h ||_{L^2(0,1)}', ...
               'L2 error vs dof (adaptive)', false);
end
if prm.plots.ieff
  plot_metric (R, prm, 'uniform',  'Ieff', 'I_{eff} = \eta / || u - u_h ||_{L^2}', ...
               'I_eff vs dof (uniform)', true);
  plot_metric (R, prm, 'adaptive', 'Ieff', 'I_{eff} = \eta / || u - u_h ||_{L^2}', ...
               'I_eff vs dof (adaptive)', true);
end

end 

%%  THE ONE REFINEMENT LOOP  (uniform and adaptive share it)

function S = run_1d (p, reg, nsteps, mode, prm)
%RUN_1D  Solve -> estimate -> mark -> insert knots.
geometry = geo_load (nrbline ([0 0], [1 0]));
zeta = linspace (0, 1, prm.nsub0 + 1);

S.dof = [];  S.L2 = [];  S.eta = [];  S.Ieff = [];

for it = 1:nsteps
  [msh_h, sp_h, knots] = build_space (geometry, zeta, p, reg, p + prm.nquad_add);
  if prm.twin_quad
    [msh_s, sp_s] = build_space (geometry, zeta, p, reg, p + prm.nquad_stab_add);
  else
    msh_s = msh_h;  sp_s = sp_h;
  end

  [K, F] = assemble (sp_s, msh_s, sp_h, msh_h, zeta, prm);
  u = solve_dirichlet (K, F, sp_h);

  E = estimate_and_error (knots, p, zeta, u, prm);

  S.dof(end+1)  = sp_h.ndof;                                  %#ok<AGROW>
  S.L2(end+1)   = E.err;                                      %#ok<AGROW>
  S.eta(end+1)  = E.eta;                                      %#ok<AGROW>
  S.Ieff(end+1) = E.eta / max (E.err, realmin);                %#ok<AGROW>

  if it == nsteps, break; end
  if sp_h.ndof >= prm.max_dof, break; end
  if E.eta <= prm.tol, break; end

  if strcmpi (mode, 'uniform')
    sel = 1:(numel (zeta) - 1);
  else
    sel = mark_elements (E.eta_K, prm);
  end
  if isempty (sel), break; end
  zeta = insert_knots (zeta, sel);
end
end

function zeta = insert_knots (zeta, sel)
%INSERT_KNOTS  Bisect selected knot spans.
mid  = 0.5 * (zeta(sel) + zeta(sel+1));
zeta = unique ([zeta(:).', mid(:).']);
end

function sel = mark_elements (est, prm)
est = est(:);
if all (est <= 0), sel = []; return; end
switch lower (prm.mark_strategy)
  case 'maximum'
    sel = find (est >= prm.mark_param * max (est));
  case 'dorfler'
    [sv, ord] = sort (est, 'descend');
    cs = cumsum (sv.^2);
    nt = find (cs >= prm.mark_param * cs(end), 1, 'first');
    if isempty (nt), nt = numel (sv); end
    sel = ord(1:nt);
  otherwise
    error ('iga:mark', 'Unknown mark_strategy ''%s''.', prm.mark_strategy);
end
sel = sort (sel(:).');
end


%% ======================================================================
%  SPACE, ASSEMBLY, SOLVE
%  ======================================================================
function [msh, space, knots] = build_space (geometry, zeta, p, reg, nquad)
%BUILD_SPACE  B-spline space on breakpoints ZETA; interior multiplicity p-reg
%   (reg = p-1 gives C^{p-1}).
zeta = unique (zeta(:).');
mult = max (1, p - reg);
inner = reshape (repmat (zeta(2:end-1), mult, 1), 1, []);
kv = [repmat(zeta(1), 1, p+1), inner, repmat(zeta(end), 1, p+1)];
knots = {kv};
rule    = msh_gauss_nodes (nquad);
[qn,qw] = msh_set_quad_nodes ({zeta}, rule);
msh     = msh_cartesian ({zeta}, qn, qw, geometry);
space   = sp_bspline (knots, p, msh);
end

function [K, F] = assemble (sp_s, msh_s, sp_h, msh_h, zeta, prm)
%ASSEMBLE  Galerkin (p+5 rule) + stabilisation (msh_s rule), triplet form.
ndof = sp_h.ndof;  F = zeros (ndof, 1);
he = diff (unique (zeta(:).'));

% Galerkin: (w, a u') + (w', kappa u') + (w, s u) - (w, f)
me = msh_evaluate_element_list (msh_h, 1:msh_h.nel);
se = sp_evaluate_element_list  (sp_h, me, 'value', true, 'gradient', true);
[N, Nx, ~, ~, Jw, conn, nsh] = unpack (se, me);
[I1, J1, V1, F] = loop_galerkin (N, Nx, Jw, conn, nsh, me.nel, F, prm);

% Stabilisation: + tau ( P w , L u - f )
if prm.stabilise
  me = msh_evaluate_element_list (msh_s, 1:msh_s.nel);
  se = sp_evaluate_element_list  (sp_s, me, 'value', true, 'gradient', true, ...
                                  'hessian', true);
  [N, Nx, Nxx, ~, Jw, conn, nsh] = unpack (se, me);
  [I2, J2, V2, F] = loop_stab (N, Nx, Nxx, Jw, he, conn, nsh, me.nel, F, prm);
else
  I2 = [];  J2 = [];  V2 = [];
end

K = sparse ([I1; I2], [J1; J2], [V1; V2], ndof, ndof);
end

function [Iv, Jv, Vv, F] = loop_galerkin (N, Nx, Jw, conn, nsh, nel, F, prm)
a = prm.a; kap = prm.kappa; s = prm.s; f = prm.f;
cap = sum (nsh(1:nel).^2);
Iv = zeros (cap,1);  Jv = Iv;  Vv = Iv;  ptr = 0;
for iel = 1:nel
  m = nsh(iel);  idx = conn(1:m,iel);
  Nv = N(:,1:m,iel);  Bx = Nx(:,1:m,iel);  w = Jw(:,iel);
  Ke = (Nv'*(w.*(a*Bx))) + (Bx'*(w.*(kap*Bx))) + (Nv'*(w.*(s*Nv)));
  q = ptr + (1:m*m);
  Iv(q) = repmat (idx, m, 1);
  Jv(q) = reshape (repmat (idx.', m, 1), [], 1);
  Vv(q) = Ke(:);
  ptr = ptr + m*m;
  F(idx) = F(idx) + Nv'*(w.*f);
end
Iv = Iv(1:ptr);  Jv = Jv(1:ptr);  Vv = Vv(1:ptr);
end

function [Iv, Jv, Vv, F] = loop_stab (N, Nx, Nxx, Jw, he, conn, nsh, nel, F, prm)
%LOOP_STAB  Stabilisation: 'sgs' P=-L^*w, 'supg' P=a w', 'gls' P=Lw.
a = prm.a; kap = prm.kappa; s = prm.s; f = prm.f;
cap = sum (nsh(1:nel).^2);
Iv = zeros (cap,1);  Jv = Iv;  Vv = Iv;  ptr = 0;
for iel = 1:nel
  m = nsh(iel);  idx = conn(1:m,iel);
  Nv  = N(:,1:m,iel);  Bx = Nx(:,1:m,iel);  Bxx = Nxx(:,1:m,iel);
  w   = Jw(:,iel);   tau = tau_stab (he(iel), prm);
  Lu = a*Bx - kap*Bxx + s*Nv;
  switch lower (prm.test_op)
    case 'sgs',  P = a*Bx + kap*Bxx - s*Nv;
    case 'supg', P = a*Bx;
    case 'gls',  P = Lu;
  end
  Ke = P'*((tau.*w).*Lu);
  q = ptr + (1:m*m);
  Iv(q) = repmat (idx, m, 1);
  Jv(q) = reshape (repmat (idx.', m, 1), [], 1);
  Vv(q) = Ke(:);
  ptr = ptr + m*m;
  F(idx) = F(idx) + P'*((tau.*w).*f);
end
Iv = Iv(1:ptr);  Jv = Jv(1:ptr);  Vv = Vv(1:ptr);
end

function u = solve_dirichlet (K, F, space)
drch = [];
try
  drch = [space.boundary(1).dofs(:); space.boundary(2).dofs(:)];
catch
  drch = [];
end
if isempty (drch), drch = [1; space.ndof]; end     % 1D open knot vector
drch = unique (drch);
free = setdiff (1:space.ndof, drch);
u = zeros (space.ndof, 1);
u(free) = K(free,free) \ F(free);
end


%%  tau

function tau = tau_stab (he, prm)
% tau_stab = (h/2|a|)(coth Pe - 1/Pe),  Pe = |a|h/(2 kappa)
a  = abs (prm.a);
Pe = a * he / (2 * prm.kappa);
tau = he / (2 * a) * xi_fun (Pe);
end

function tau = tau_est (he, prm)
% tau_est = min( h/(sqrt(2)|a|) , h^2/(3 sqrt(10) kappa) )
a = abs (prm.a);
tau = min (he / (sqrt (2) * a), he^2 / (3 * sqrt (10) * prm.kappa));
end

function xi = xi_fun (Pe)
if Pe < 1e-8, xi = Pe / 3; else, xi = coth (Pe) - 1 / Pe; end
end


%%  ESTIMATOR + TRUE ERROR  

function E = estimate_and_error (knots, p, zeta, u, prm)
%ESTIMATE_AND_ERROR  Per-element eta_K (VMS residual) and ||e||_K.
%   Elements near the boundary layer are split into graded sub-cells.
a = prm.a; kap = prm.kappa; s = prm.s; f = prm.f;
zeta = unique (zeta(:).');
nel  = numel (zeta) - 1;
sub  = strcmpi (prm.quad_mode, 'subcell');
band = prm.lay_bandw * prm.delta;
nq   = p + prm.nquad_add;

etaK = zeros (nel,1);  errK = zeros (nel,1);

for iel = 1:nel
  aa = zeta(iel);  bb = zeta(iel+1);  h = bb - aa;
  ns = 1;  gr = [];
  if sub && bb > 1 - band && h > prm.delta/2
    ns = min (max (ceil (h / (0.25*prm.delta)), 1), prm.lay_nsub);
    gr = prm.lay_grade;
  end
  [x, w] = subcell_rule (aa, bb, nq, ns, gr);
  B = colloc_ders (knots{1}, p, x, 2);
  uh = B{1}*u;  uhx = B{2}*u;  uhxx = B{3}*u;

  errK(iel) = sqrt (max (sum (w(:) .* (prm.uex(x(:)) - uh).^2), 0));
  res       = a*uhx - kap*uhxx + s*uh - f;
  etaK(iel) = tau_est (h, prm) * sqrt (max (sum (w(:) .* res.^2), 0));
end

E.eta_K = etaK;  E.err_K = errK;
E.err = sqrt (sum (errK.^2));
E.eta = sqrt (sum (etaK.^2));
end

function [xq, wq] = subcell_rule (a, b, nq, ns, grade)
%SUBCELL_RULE  Gauss on [a,b], optionally split into NS cells graded toward b.
[xg, wg] = gauss_ref (nq);
if ns <= 1
  e = [a b];
elseif isempty (grade)
  e = linspace (a, b, ns+1);
else
  e = a + (b-a) * linspace (0, 1, ns+1).^grade;
end
xq = zeros (1, ns*nq);  wq = xq;
for i = 1:ns
  aa = e(i);  bb = e(i+1);
  q  = (i-1)*nq + (1:nq);
  xq(q) = 0.5*(aa+bb) + 0.5*(bb-aa)*xg;
  wq(q) = 0.5*(bb-aa)*wg;
end
end

function [x, w] = gauss_ref (n)
%GAUSS_REF  n-point Gauss-Legendre on [-1,1] via Golub-Welsch.
k  = 1:n-1;
bb = k ./ sqrt (4*k.^2 - 1);
T  = diag (bb,1) + diag (bb,-1);
[Vv, Dd] = eig (T);
[x, idx] = sort (diag (Dd));
w = 2 * (Vv(1,idx).^2);
x = x(:).';  w = w(:).';
end



%%  B-SPLINE EVALUATION

function B = colloc_ders (knots, p, x, nd)
%COLLOC_DERS  Collocation matrices for basis and derivatives up to order ND.
%   At p=1, derivatives above order 1 are exact zeros (so -kappa u'' = 0).
knots = knots(:).';
n = numel (knots) - p - 1;
x = min (max (x(:).', knots(1)), knots(end));
nx = numel (x);
sp = findspan (n-1, p, x, knots);
nu = min (nd, p);
D  = basisfunder (sp, p, x, knots, nu);
rows = repmat ((1:nx).', 1, p+1);
cols = bsxfun (@plus, double (sp(:)) - p, 0:p) + 1;
cols = min (max (cols, 1), n);
B = cell (1, nd+1);
for d = 0:nd
  if d > nu
    B{d+1} = sparse (nx, n);
  else
    Vd = reshape (D(:, d+1, :), nx, p+1);
    B{d+1} = sparse (rows(:), cols(:), Vd(:), nx, n);
  end
end
end

function [N, Nx, Nxx, xq, Jw, conn, nsh, hsz] = unpack (sp, msh)
nqn = msh.nqn;  nel = msh.nel;  nshmax = sp.nsh_max;
N = reshape (sp.shape_functions, nqn, nshmax, nel);
if isfield (sp, 'shape_function_gradients')
  G  = sp.shape_function_gradients;
  Nx = reshape (G(1,:,:,:), nqn, nshmax, nel);
else
  Nx = zeros (nqn, nshmax, nel);
end
if isfield (sp, 'shape_function_hessians')
  H   = sp.shape_function_hessians;
  Nxx = reshape (H(1,1,:,:,:), nqn, nshmax, nel);
else
  Nxx = zeros (nqn, nshmax, nel);
end
xq   = reshape (msh.geo_map(1,:,:), nqn, nel);
Jw   = msh.jacdet .* msh.quad_weights;
conn = sp.connectivity;
nsh  = sp.nsh;
hsz  = (sum (Jw, 1)).^(1/msh.ndim);
end



%%  SUMMARY

function print_summary (R, prm)
ok = find ([R.ok]);
if isempty (ok), return; end
fprintf ('\n%s\n', repmat ('=', 1, 76));
fprintf (' SUMMARY   (kappa = %g, test_op = ''%s'')\n', prm.kappa, lower (prm.test_op));
fprintf ('%s\n', repmat ('=', 1, 76));
fprintf ('  %2s  %-9s  %8s  %10s  %10s  %8s\n', ...
         'p', 'mode', 'dof', 'eta', 'L2', 'Ieff');
fprintf ('  %s\n', repmat ('-', 1, 72));
for k = ok
  fprintf ('  %2d  %-9s  %8d  %10.4e  %10.4e  %8.3f\n', ...
           R(k).degree, R(k).mode, R(k).dof(end), R(k).eta(end), ...
           R(k).L2(end), R(k).Ieff(end));
end
fprintf ('%s\n\n', repmat ('=', 1, 76));
end



%%  PLOTS

function C = degree_colors ()
C = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
     0.49 0.18 0.56; 0.93 0.69 0.13];
end

function M = degree_markers ()
M = {'o','s','^','d','v'};
end

function plot_metric (R, prm, mode, field, ylab, figname, is_ieff)
%PLOT_METRIC  FIELD vs dof for every degree, restricted to MODE.
sel = false (1, numel (R));
for k = 1:numel (R)
  sel(k) = R(k).ok && strcmpi (R(k).mode, mode);
end
if ~any (sel), return; end
figure ('Name', figname, 'Color','w');
C = degree_colors ();  M = degree_markers ();
hold on;  leg = {};
for k = find (sel)
  p = R(k).degree;
  plot (R(k).dof, R(k).(field), [M{p} '-'], 'Color', C(p,:), ...
        'LineWidth',1.7,'MarkerFaceColor',C(p,:),'MarkerSize',6);
  leg{end+1} = sprintf ('p = %d, C^{%d}', p, R(k).regularity);   %#ok<AGROW>
end
set (gca,'XScale','log','YScale','log');
if nargin > 6 && is_ieff
  xl = xlim;  yl = ylim;
  patch ([xl(1) xl(2) xl(2) xl(1)], [1 1 10 10], [0.85 0.92 0.85], ...
         'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
  plot (xl, [1 1], 'k-', 'LineWidth',1.3,'HandleVisibility','off');
  try, uistack (findobj (gca,'Type','patch'), 'bottom'); catch, end
  xlim (xl);  ylim (yl);
end
grid on; box on;
xlabel ('degrees of freedom');  ylabel (ylab);
if ~isempty (leg), legend (leg,'Location','southwest','FontSize',8); end
hold off;
end


%%  OPTION PARSING

function prm = parse_options (prm, args)
if numel (args) == 1 && isstruct (args{1})
  f = fieldnames (args{1});
  for i = 1:numel (f), prm.(f{i}) = args{1}.(f{i}); end
  return;
end
if mod (numel (args), 2) ~= 0
  error ('iga:args', 'Options must come in name/value pairs.');
end
for i = 1:2:numel (args)
  name = args{i};
  if ~ischar (name), error ('iga:args', 'Option names must be strings.'); end
  if any (strcmpi (name, {'degree','regularity'}))
    error ('iga:args', ...
      ['''%s'' no longer exists: the five pairs (1,C^0)...(5,C^4) are fixed. ' ...
       'Use ''degrees'', e.g. ''degrees'',[1 3].'], name);
  end
  if any (strcmpi (name, {'estimator','estimators','drivers','gradind_scale'}))
    error ('iga:args', ...
      ['''%s'' no longer exists: the gradient indicator has been removed and ' ...
       'the VMS estimator is the only one computed.'], name);
  end
  if any (strcmpi (name, {'refinement'}))
    error ('iga:args', ...
      ['''%s'' no longer exists: uniform and adaptive refinement both run ' ...
       'automatically every call.'], name);
  end
  if isfield (prm, name)
    prm.(name) = args{i+1};
  elseif isfield (prm.plots, name)
    prm.plots.(name) = args{i+1};
  else
    error ('iga:args', 'Unknown option ''%s''.', name);
  end
end
end