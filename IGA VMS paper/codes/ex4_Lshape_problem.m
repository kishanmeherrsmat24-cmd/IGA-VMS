function out = iga_thb_aniso_Lshape (varargin)
%IGA_THB_ANISO_LSHAPE  Anisotropic THB-splines on 3-patch L-shaped domain.
%     -eps Lap(u)+b.grad(u)+sigma u=f in Omega,  u=0 on dOmega
%     Omega=(0,1)^2\[0.5,1]x[0,0.5], eps=1e-6, b=(1,3)^T, sigma=1
%   Two drivers: 'vms' (eta_K=tau_K||R||) and 'gradind' (||grad u_h||_K).
%   Four knot vectors (interfaces conforming at every level):
%     patch 1=[0,0.5]x[0,0.5]: ku=U12, kv=V1
%     patch 2=[0,0.5]x[0.5,1]: ku=U12, kv=V23
%     patch 3=[0.5,1]x[0.5,1]: ku=U3,  kv=V23
%   No exact solution, so no L2 error or I_eff reported.
%   Requires GeoPDEs 3.x (multipatch) and the NURBS toolbox.

%%                          USER PARAMETERS
prm.eps       = 1e-6;
prm.b         = [1; 3];
prm.c         = 1;
prm.degrees   = [1 2 3 4 5];     % -> p = 1..5, highest continuity C^{p-1}
prm.drivers   = {'vms','gradind'};
prm.nel0      = 7;
prm.niter     = 40;
prm.max_ndof  = 4000;
prm.theta     = 0.15;
prm.mark_power= 1;
prm.dir_kappa = 0.50;
prm.isotropic = false;
prm.test_op   = 'sgs';
prm.stabilise = true;
prm.norm_type = 'hs';
prm.gradind_scale = 'none';      % 'none' | 'h'
prm.nquad_add = 2;
prm.gates     = true;
prm.verbose   = true;
% --- figures ----------------------------------------------------------
prm.plot_npts = 90;              % points per direction per patch, solution
prm.plots.eta      = true;       % eta vs dof, all runs, one axes
prm.plots.mesh     = true;       % final mesh, ONE FIGURE PER (degree,driver)
prm.plots.solution = true;       % solution,  ONE FIGURE PER (degree,driver)
prm.fig_export= '';              % folder for pdf/png output; '' = do not save
prm.fig_fmt   = {'pdf','png'};
prm.fig_font  = 'Times New Roman';
prm.fig_fsize = 11;
prm.xpad      = 2.2;
% ======================================================================
prm = parse_options (prm, varargin);
if ischar (prm.drivers), prm.drivers = {prm.drivers}; end
if ischar (prm.fig_fmt), prm.fig_fmt = {prm.fig_fmt}; end

R = struct ('degree',{},'driver',{},'hist',{},'H',{},'u',{},'ok',{});
m = 0;
for j = 1:numel (prm.drivers)
  for k = 1:numel (prm.degrees)
    m = m + 1;
    R(m).degree = prm.degrees(k);  R(m).driver = prm.drivers{j};  R(m).ok = false;
    try
      [hist, H, u] = run_one (prm.degrees(k), prm.drivers{j}, prm);
      R(m).hist = hist;  R(m).H = H;  R(m).u = u;  R(m).ok = true;
    catch err
      fprintf (2, '  p = %d, driver ''%s'' FAILED: %s\n', prm.degrees(k), ...
               prm.drivers{j}, err.message);
      if ~isempty (err.stack)
        for q = 1:min (3, numel (err.stack))
          fprintf (2, '     at %s line %d\n', err.stack(q).name, err.stack(q).line);
        end
      end
    end
  end
end
out = struct ('runs', {R}, 'prm', prm);
print_tables (R, prm);
if prm.plots.eta,      plot_eta_all      (R, prm); end
if prm.plots.mesh,     plot_final_meshes (R, prm); end
if prm.plots.solution, plot_solution_surfs (R, prm); end
if nargout == 0, clear out; end
end 


function [hist, H, u] = run_one (p, driver, prm)
pb = lshape_problem (prm);
P  = lshape_patches ();
t  = linspace (0,1,prm.nel0+1);
k0 = open_knots (t, p);
H.p  = p;
H.P  = P;
H.KN = { struct('U12',k0,'V1',k0,'V23',k0,'U3',k0) };
H.Om = { [0 1 0 1] };
H = build_level (H, 1, prm);
H = rebuild_sets (H);
H = update_functions (H);
H = build_Csub (H);

hist = struct ('ndof',[],'nel',[],'eta',[],'eta_vms',[],'eta_grad',[], ...
               'nlev',[],'aspect',[],'tsec',[]);
u = [];
for it = 1:prm.niter
  t0 = tic;
  if prm.gates, check_gates (H, it); end
  u = solve_stab (H, pb, prm);
  E = estimate (H, u, pb, prm);
  if strcmpi (driver,'gradind'), ed = E.eta_grad; else, ed = E.eta; end
  hist.ndof(end+1)     = H.ndof;
  hist.nel(end+1)      = numel (E.eta_K);
  hist.eta(end+1)      = ed;
  hist.eta_vms(end+1)  = E.eta;
  hist.eta_grad(end+1) = E.eta_grad;
  hist.nlev(end+1)     = numel (H.KN);
  hist.aspect(end+1)   = E.aspect;
  hist.tsec(end+1)     = toc (t0);
  if prm.verbose
    fprintf (['  p=%d %-8s it %2d  lev %2d  ndof %6d  nel %6d  eta %.4e' ...
              '  aspect %.1e  %.1fs\n'], p, driver, it, numel (H.KN), ...
             H.ndof, hist.nel(end), ed, E.aspect, hist.tsec(end));
  end
  if it == prm.niter || H.ndof >= prm.max_ndof, break; end
  H = refine (H, E, driver, prm);
end
end


%%  PROBLEM AND PATCHES

function pb = lshape_problem (prm)
pb.eps = prm.eps;  pb.b = prm.b(:);  pb.c = prm.c;
pb.bnorm = norm (pb.b);  pb.ahat = pb.b/pb.bnorm;
pb.f = @(x,y) lshape_rhs (x,y);
end
function v = lshape_rhs (x, y)
r = sqrt ((x-0.5).^2 + (y-0.5).^2);
v = 100 * r .* (r - 0.5) .* (r - sqrt (0.5));
end

function P = lshape_patches ()
P(1) = struct ('x0',0.0,'Lx',0.5,'y0',0.0,'Ly',0.5,'ku','U12','kv','V1', ...
               'ext',[1 2 3]);
P(2) = struct ('x0',0.0,'Lx',0.5,'y0',0.5,'Ly',0.5,'ku','U12','kv','V23', ...
               'ext',[1 4]);
P(3) = struct ('x0',0.5,'Lx',0.5,'y0',0.5,'Ly',0.5,'ku','U3', 'kv','V23', ...
               'ext',[2 3 4]);
end

function [bnd, itf] = lshape_topology ()
itf(1).patch1 = 1;  itf(1).side1 = 4;
itf(1).patch2 = 2;  itf(1).side2 = 3;  itf(1).ornt = 1;
itf(2).patch1 = 2;  itf(2).side1 = 2;
itf(2).patch2 = 3;  itf(2).side2 = 1;  itf(2).ornt = 1;
sp = {1,1; 1,2; 1,3; 2,1; 2,4; 3,2; 3,3; 3,4};
for k = 1:size (sp,1)
  bnd(k).nsides = 1;  bnd(k).patches = sp{k,1};  bnd(k).faces = sp{k,2}; %#ok<AGROW>
end
end

function k = open_knots (brk, p)
brk = unique (brk(:).');
k = [repmat(brk(1),1,p+1), brk(2:end-1), repmat(brk(end),1,p+1)];
end


%%  LEVEL CONSTRUCTION

function H = build_level (H, l, prm)
p  = H.p;
nq = p + prm.nquad_add;
[bnd, itf] = lshape_topology ();
np = numel (H.P);
msh_a = cell (1,np);  sp_a = cell (1,np);
for q = 1:np
  Pq = H.P(q);
  kn = {H.KN{l}.(Pq.ku), H.KN{l}.(Pq.kv)};
  brk = {unique(kn{1}), unique(kn{2})};
  c = zeros (2,2,2);
  x0 = Pq.x0; x1 = Pq.x0+Pq.Lx; y0 = Pq.y0; y1 = Pq.y0+Pq.Ly;
  c(1,1,1)=x0; c(1,2,1)=x1; c(1,1,2)=x0; c(1,2,2)=x1;
  c(2,1,1)=y0; c(2,2,1)=y0; c(2,1,2)=y1; c(2,2,2)=y1;
  geo = geo_load (nrbmak (c, {[0 0 1 1],[0 0 1 1]}));
  [qn, qw] = msh_set_quad_nodes (brk, msh_gauss_nodes ([nq nq]));
  msh_a{q} = msh_cartesian (brk, qn, qw, geo);
  sp_a{q}  = sp_bspline (kn, [p p], msh_a{q});
  H.brk{l}{q} = brk;
  H.nd{l}{q}  = [numel(brk{1})-1, numel(brk{2})-1];
end
msh_mp = msh_multipatch (msh_a, bnd);
H.msh{l}  = msh_a;
H.sp{l}   = sp_a;
H.spmp{l} = sp_multipatch (sp_a, msh_mp, itf);
H.nmp(l)  = H.spmp{l}.ndof;
end

function tf = in_domain (Om, cx, cy)
tf = false (numel (cx), 1);
if isempty (Om), return; end
tol = 1e-12;
for k = 1:size (Om,1)
  tf = tf | (cx >= Om(k,1)-tol & cx <= Om(k,2)+tol & ...
             cy >= Om(k,3)-tol & cy <= Om(k,4)+tol);
end
end

function [cx, cy] = cell_centres (H, l, q)
Pq = H.P(q);  brk = H.brk{l}{q};
bx = Pq.x0 + Pq.Lx * brk{1}(:);
by = Pq.y0 + Pq.Ly * brk{2}(:);
cxv = 0.5*(bx(1:end-1)+bx(2:end));
cyv = 0.5*(by(1:end-1)+by(2:end));
[CX, CY] = ndgrid (cxv, cyv);
cx = CX(:);  cy = CY(:);
end

function H = rebuild_sets (H)
L = numel (H.KN);
for l = 1:L
  for q = 1:numel (H.P)
    [cx, cy] = cell_centres (H, l, q);
    inl = in_domain (H.Om{l}, cx, cy);
    if l < L, nxt = in_domain (H.Om{l+1}, cx, cy); else, nxt = false (size (inl)); end
    H.act{l}{q} = find (inl & ~nxt);
    H.dea{l}{q} = find (inl &  nxt);
  end
end
end

function H = update_functions (H)
%UPDATE_FUNCTIONS  Active/deactivated multipatch functions; interface functions use AND/OR.
L = numel (H.KN);  p = H.p;
for l = 1:L
  n = H.nmp(l);
  allOm = true (n,1);  anyAct = false (n,1);  seen = false (n,1);
  for q = 1:numel (H.P)
    nd = H.nd{l}{q};
    n1 = numel (H.sp{l}{q}.knots{1}) - p - 1;
    n2 = numel (H.sp{l}{q}.knots{2}) - p - 1;
    inOm = false (nd);  inOm([H.act{l}{q}(:); H.dea{l}{q}(:)]) = true;
    inAc = false (nd);  inAc(H.act{l}{q}(:)) = true;
    SO = cumsum (cumsum (double (inOm),1),2);
    SA = cumsum (cumsum (double (inAc),1),2);
    g  = H.spmp{l}.gnum{q}(:);
    for j = 1:n2
      j0 = max (1, j-p);  j1 = min (nd(2), j);
      for i = 1:n1
        i0 = max (1, i-p);  i1 = min (nd(1), i);
        f  = i + (j-1)*n1;
        F  = g(f);
        ncell = (i1-i0+1)*(j1-j0+1);
        allOm(F)  = allOm(F)  && (block_sum (SO,i0,i1,j0,j1) >= ncell-0.5);
        anyAct(F) = anyAct(F) || (block_sum (SA,i0,i1,j0,j1) > 0.5);
        seen(F)   = true;
      end
    end
  end
  H.actf{l} = find (seen & allOm &  anyAct);
  H.deaf{l} = find (seen & allOm & ~anyAct);
end
end

function H = build_Csub (H)
L = numel (H.KN);  p = H.p;
C = cell (1,L);
C{1} = sparse (H.actf{1}, 1:numel (H.actf{1}), 1, H.nmp(1), numel (H.actf{1}));
for l = 2:L
  I = [];  J = [];  V = [];
  for q = 1:numel (H.P)
    Pq = H.P(q);
    Px = basiskntins (p, H.KN{l-1}.(Pq.ku), H.KN{l}.(Pq.ku));
    Py = basiskntins (p, H.KN{l-1}.(Pq.kv), H.KN{l}.(Pq.kv));
    Tq = kron (sparse (Py), sparse (Px));
    [a, b, v] = find (Tq);
    gf = H.spmp{l}.gnum{q}(:);      gc = H.spmp{l-1}.gnum{q}(:);
    I = [I; gf(a)];  J = [J; gc(b)];  V = [V; v];   %#ok<AGROW>
  end
  % interface functions are written twice with the SAME value: de-duplicate
  % rather than let sparse() add them together.
  [IJ, ia] = unique ([I J], 'rows');
  T = sparse (IJ(:,1), IJ(:,2), V(ia), H.nmp(l), H.nmp(l-1));
  kill = union (H.actf{l}, H.deaf{l});
  T(kill,:) = 0;                                     % TRUNCATION
  Ecol = sparse (H.actf{l}, 1:numel (H.actf{l}), 1, H.nmp(l), numel (H.actf{l}));
  C{l} = [T * C{l-1}, Ecol];
end
H.C = C;
H.ndof = size (C{L}, 2);
end

function s = block_sum (S, i0, i1, j0, j1)
s = S(i1,j1);
if i0 > 1, s = s - S(i0-1,j1); end
if j0 > 1, s = s - S(i1,j0-1); end
if i0 > 1 && j0 > 1, s = s + S(i0-1,j0-1); end
end

function check_gates (H, it)
L = numel (H.KN);
nm = {'U12','V1','V23','U3'};
for l = 2:L
  for d = 1:4
    a = round (H.KN{l-1}.(nm{d})*1e12)/1e12;
    b = round (H.KN{l}.(nm{d})*1e12)/1e12;
    if ~all (ismember (a, b))
      error ('thb:nesting','GATE FAIL (it %d): %s not nested, level %d->%d.', ...
             it, nm{d}, l-1, l);
    end
  end
end
for l = 1:L
  if all (cellfun (@isempty, H.act{l})), continue; end
  s = H.C{l} * ones (size (H.C{l},2), 1);
  for q = 1:numel (H.P)
    if isempty (H.act{l}{q}), continue; end
    g = H.spmp{l}.gnum{q}(:);
    p = H.p;  nd = H.nd{l}{q};
    n1 = numel (H.sp{l}{q}.knots{1}) - p - 1;
    n2 = numel (H.sp{l}{q}.knots{2}) - p - 1;
    [ia, ja] = ind2sub (nd, H.act{l}{q}(:));
    % cell (i,j) lies in the support of the tensor functions i:i+p, j:j+p.
    % Partition of unity means each of their COEFFICIENTS is 1, since the
    % tensor functions themselves already sum to 1 - not that they sum to 1.
    fl = [];
    for dx = 0:p
      for dy = 0:p
        fl = [fl; min(ia+dx,n1) + (min(ja+dy,n2)-1)*n1];   %#ok<AGROW>
      end
    end
    F = unique (g(fl));
    e = max (abs (s(F) - 1));
    if e > 1e-8
      error ('thb:pou',['GATE FAIL (it %d): partition of unity off by ' ...
             '%.2e at level %d, patch %d.'], it, e, l, q);
    end
  end
end
end

%%  SOLVER  (SGS stabilisation )

function u = solve_stab (H, pb, prm)
N = H.ndof;
A = sparse (N,N);  F = zeros (N,1);
for l = 1:numel (H.KN)
  nl = H.nmp(l);  Nl = size (H.C{l},2);
  Al = sparse (nl,nl);  Fl = zeros (nl,1);
  any_el = false;
  for q = 1:numel (H.P)
    if isempty (H.act{l}{q}), continue; end
    any_el = true;
    [msh_el, sp_el] = eval_cells (H, l, q, H.act{l}{q}, false);
    [hx, hy] = elem_sizes (H, l, q, H.act{l}{q});
    [Ke, Fe, cn] = elem_matrices (msh_el, sp_el, hx, hy, pb, prm);
    g  = H.spmp{l}.gnum{q}(:);
    cn = g(cn);
    ns = size (cn,1);  ne = size (cn,2);
    I = reshape (repmat (reshape (cn,[ns 1 ne]), [1 ns 1]), [], 1);
    J = reshape (repmat (reshape (cn,[1 ns ne]), [ns 1 1]), [], 1);
    Al = Al + sparse (I, J, Ke(:), nl, nl);
    Fl = Fl + accumarray (cn(:), reshape (Fe,[],1), [nl 1]);
  end
  if ~any_el, continue; end
  A(1:Nl,1:Nl) = A(1:Nl,1:Nl) + H.C{l}.' * Al * H.C{l};
  F(1:Nl)      = F(1:Nl)      + H.C{l}.' * Fl;
end
drch = boundary_dofs (H);
free = setdiff ((1:N).', drch(:));
Af = A(free,free);  bf = F(free);
dv = full (sqrt (abs (diag (Af))));
dv(~isfinite (dv) | dv <= 0) = 1;
D = spdiags (1./dv, 0, numel (dv), numel (dv));
u = zeros (N,1);
u(free) = D * ((D*Af*D) \ (D*bf));
end

function d = boundary_dofs (H)
d = [];  off = 0;
for l = 1:numel (H.KN)
  bd = [];
  for q = 1:numel (H.P)
    g = H.spmp{l}.gnum{q}(:);
    for s = H.P(q).ext
      bd = union (bd, g(H.sp{l}{q}.boundary(s).dofs(:)));
    end
  end
  tf = ismember (H.actf{l}, bd);
  ff = find (tf);
  d = [d; off + ff];                 %#ok<AGROW>
  off = off + numel (H.actf{l});
end
d = unique (d);
end

function [msh_el, sp_el] = eval_cells (H, l, q, elems, want_hess)
msh_el = msh_evaluate_element_list (H.msh{l}{q}, elems);
if want_hess
  sp_el = sp_evaluate_element_list (H.sp{l}{q}, msh_el, 'value',true, ...
              'gradient',true,'hessian',true);
  Hs = sp_el.shape_function_hessians;
  sp_el.shape_function_laplacians = ...
      reshape (Hs(1,1,:,:,:)+Hs(2,2,:,:,:), size(Hs,3), size(Hs,4), size(Hs,5));
else
  sp_el = sp_evaluate_element_list (H.sp{l}{q}, msh_el, 'value',true, ...
              'gradient',true,'laplacian',true);
end
end

function [hx, hy] = elem_sizes (H, l, q, elems)
Pq = H.P(q);  brk = H.brk{l}{q};
bx = Pq.x0 + Pq.Lx*brk{1}(:);
by = Pq.y0 + Pq.Ly*brk{2}(:);
[i1, j1] = ind2sub (H.nd{l}{q}, double (elems(:)));
hx = bx(i1+1)-bx(i1);   hy = by(j1+1)-by(j1);
end

function [Ke, Fe, conn] = elem_matrices (msh_el, sp_el, hx, hy, pb, prm)
[nqn, nsh, nel] = size (sp_el.shape_functions);
conn = double (sp_el.connectivity);
W  = msh_el.quad_weights .* msh_el.jacdet;
W3 = reshape (W, nqn, 1, nel);
N  = sp_el.shape_functions;
Gx = reshape (sp_el.shape_function_gradients(1,:,:,:), nqn, nsh, nel);
Gy = reshape (sp_el.shape_function_gradients(2,:,:,:), nqn, nsh, nel);
Lp = sp_el.shape_function_laplacians;
x  = reshape (msh_el.geo_map(1,:,:), nqn, nel);
y  = reshape (msh_el.geo_map(2,:,:), nqn, nel);
fq = pb.f (x,y);
bG  = pb.b(1)*Gx + pb.b(2)*Gy;
Rop = -pb.eps*Lp + bG + pb.c*N;
if prm.stabilise
  switch lower (prm.test_op)
    case 'sgs',  P = bG + pb.eps*Lp - pb.c*N;
    case 'supg', P = bG;
    otherwise,   P = Rop;
  end
else
  P = zeros (nqn,nsh,nel);
end
tau3 = reshape (tau_stab (hx, hy, pb), 1, 1, nel);
NT = permute (N,[2 1 3]);
Ke = pagemtimes (NT, W3.*bG) ...
   + pb.eps*(pagemtimes (permute (Gx,[2 1 3]), W3.*Gx) ...
           + pagemtimes (permute (Gy,[2 1 3]), W3.*Gy)) ...
   + pb.c*pagemtimes (NT, W3.*N) ...
   + tau3 .* pagemtimes (permute (P,[2 1 3]), W3.*Rop);
Wf = reshape (W.*fq, nqn, 1, nel);
Fe = pagemtimes (NT, Wf) + tau3 .* pagemtimes (permute (P,[2 1 3]), Wf);
end

function tau = tau_stab (hx, hy, pb)
a1 = abs (pb.ahat(1));  a2 = abs (pb.ahat(2));
hf = min (hx/max(a1,1e-14), hy/max(a2,1e-14));
al = pb.bnorm * hf / (2*pb.eps);
xi = zeros (size (al));
sm = al < 1e-8;  xi(sm) = al(sm)/3;
lg = al > 50;    xi(lg) = 1 - 1./al(lg);
md = ~sm & ~lg;  xi(md) = coth (al(md)) - 1./al(md);
tau = hf ./ (2*pb.bnorm) .* xi;
end



%%  ESTIMATOR (4.13) AND GRADIENT INDICATOR (4.16)

function E = estimate (H, u, pb, prm)
%ESTIMATE  Both indicators on the same quadrature:
%   E.eta_K = tau_K ||R||_{L2(K)} (4.13),  E.eta_gradK = ||grad u_h||_{L2(K)} (4.16).
eta = [];  egr = [];  ax = [];  ay = [];  gx1 = [];  gy1 = [];
lev = [];  cel = [];  ptc = [];  hxa = [];  hya = [];  cxa = [];  cya = [];
for l = 1:numel (H.KN)
  Nl = size (H.C{l},2);
  uTP = H.C{l} * u(1:Nl);
  for q = 1:numel (H.P)
    if isempty (H.act{l}{q}), continue; end
    g = H.spmp{l}.gnum{q}(:);
    [msh_el, sp_el] = eval_cells (H, l, q, H.act{l}{q}, true);
    Hs = sp_el.shape_function_hessians;
    [nqn, nsh, nel] = size (sp_el.shape_functions);
    cn = double (sp_el.connectivity);
    ue = reshape (uTP(g(cn)), nsh, 1, nel);
    W  = msh_el.quad_weights .* msh_el.jacdet;
    W3 = reshape (W, nqn, 1, nel);
    Nb = sp_el.shape_functions;
    Gx = reshape (sp_el.shape_function_gradients(1,:,:,:), nqn, nsh, nel);
    Gy = reshape (sp_el.shape_function_gradients(2,:,:,:), nqn, nsh, nel);
    Hxx = reshape (Hs(1,1,:,:,:), nqn, nsh, nel);
    Hyy = reshape (Hs(2,2,:,:,:), nqn, nsh, nel);
    xq = reshape (msh_el.geo_map(1,:,:), nqn, nel);
    yq = reshape (msh_el.geo_map(2,:,:), nqn, nel);
    uh  = pagemtimes (Nb, ue);
    gx  = pagemtimes (Gx, ue);
    gy  = pagemtimes (Gy, ue);
    uxx = pagemtimes (Hxx, ue);
    uyy = pagemtimes (Hyy, ue);
    Rk = -pb.eps*(uxx+uyy) + pb.b(1)*gx + pb.b(2)*gy + pb.c*uh ...
         - reshape (pb.f (xq,yq), nqn, 1, nel);
    r2 = reshape (sum (W3 .* Rk.^2,  1), [], 1);
    g2 = reshape (sum (W3 .* (gx.^2 + gy.^2), 1), [], 1);
    x2 = reshape (sum (W3 .* uxx.^2, 1), [], 1);
    y2 = reshape (sum (W3 .* uyy.^2, 1), [], 1);
    px2 = reshape (sum (W3 .* gx.^2, 1), [], 1);
    py2 = reshape (sum (W3 .* gy.^2, 1), [], 1);
    [hx, hy] = elem_sizes (H, l, q, H.act{l}{q});
    [cxv, cyv] = cell_centres (H, l, q);
    tauK = tau_up (hx, hy, pb, prm);
    e_c = tauK .* sqrt (max (r2,0));
    if strcmpi (prm.gradind_scale, 'h'), sc = sqrt (hx.*hy); else, sc = ones (size (hx)); end
    g_c = sc .* sqrt (max (g2,0));
    a_c = hx.^2 .* sqrt (max (x2,0));
    b_c = hy.^2 .* sqrt (max (y2,0));
    a_g = hx    .* sqrt (max (px2,0));
    b_g = hy    .* sqrt (max (py2,0));
    eta = [eta; e_c];   egr = [egr; g_c];                       %#ok<AGROW>
    ax  = [ax; a_c];    ay  = [ay; b_c];                        %#ok<AGROW>
    gx1 = [gx1; a_g];   gy1 = [gy1; b_g];                       %#ok<AGROW>
    lev = [lev; l*ones(nel,1)];  cel = [cel; H.act{l}{q}(:)];   %#ok<AGROW>
    ptc = [ptc; q*ones(nel,1)];                                 %#ok<AGROW>
    hxa = [hxa; hx];  hya = [hya; hy];                          %#ok<AGROW>
    cxa = [cxa; cxv(H.act{l}{q})];  cya = [cya; cyv(H.act{l}{q})]; %#ok<AGROW>
  end
end
E.eta_K = eta;  E.eta_gradK = egr;
E.ax = ax;  E.ay = ay;  E.gx1 = gx1;  E.gy1 = gy1;
E.lev = lev;  E.cel = cel;  E.patch = ptc;
E.hx = hxa;  E.hy = hya;  E.cx = cxa;  E.cy = cya;
E.eta      = sqrt (sum (eta.^2));
E.eta_grad = sqrt (sum (egr.^2));
E.aspect = max ([hxa; hya]) / max (min ([hxa; hya]), realmin);
end

function t = tau_up (hx, hy, pb, prm)
switch lower (prm.norm_type)
  case 'hs', CA = 1/sqrt(2);  CD = 1/(3*sqrt(10));
  otherwise, CA = 2/pi;       CD = 1/pi^2;
end
a1 = abs (pb.ahat(1));  a2 = abs (pb.ahat(2));
ell  = 2 ./ (a1./hx + a2./hy);
ellm = min (hx/max(a1,1e-14), hy/max(a2,1e-14));
t = min (max (CA*ell/pb.bnorm, realmin), max (CD*ellm.^2/pb.eps, realmin));
end


%%  ANISOTROPIC HIERARCHICAL REFINEMENT

function H = refine (H, E, driver, prm)
%REFINE  Doerfler on functions; direction from aggregated Hessian weights;
%   inserts midpoints into the owning patch's knot vector for that direction.
if strcmpi (driver, 'gradind'), indK = E.eta_gradK; else, indK = E.eta_K; end
L = numel (H.KN);  p = H.p;
fe = [];  fl = [];  fq = [];  fi = [];  fj = [];  fx = [];  fy = [];
for l = 1:L
  for q = 1:numel (H.P)
    if isempty (H.act{l}{q}), continue; end
    nd = H.nd{l}{q};
    eg = zeros (nd);  ag = zeros (nd);  bg = zeros (nd);
    m = (E.lev == l) & (E.patch == q);
    eg(E.cel(m)) = indK(m);
    % Hessian weights, with the first-derivative weights substituted where
    % the Hessian ones vanish (p = 1: a bilinear has u_xx = u_yy = 0).
    aa = E.ax(m);  bb = E.ay(m);
    flat = (aa + bb) <= 0;
    if any (flat)
      g1 = E.gx1(m);  g2 = E.gy1(m);
      aa(flat) = g1(flat);  bb(flat) = g2(flat);
    end
    ag(E.cel(m)) = aa;
    bg(E.cel(m)) = bb;
    Se = cumsum (cumsum (eg.^prm.mark_power,1),2);
    Sa = cumsum (cumsum (ag,1),2);
    Sb = cumsum (cumsum (bg,1),2);
    n1 = numel (H.sp{l}{q}.knots{1}) - p - 1;
    n2 = numel (H.sp{l}{q}.knots{2}) - p - 1;
    g  = H.spmp{l}.gnum{q}(:);
    isact = false (H.nmp(l),1);  isact(H.actf{l}) = true;
    for j = 1:n2
      j0 = max (1,j-p);  j1 = min (nd(2), j);
      for i = 1:n1
        i0 = max (1,i-p);  i1 = min (nd(1), i);
        f = i + (j-1)*n1;
        if ~isact(g(f)), continue; end
        v = block_sum (Se,i0,i1,j0,j1);
        if v <= 0, continue; end
        fe(end+1,1) = v;   fl(end+1,1) = l;   fq(end+1,1) = q;   %#ok<AGROW>
        fi(end+1,1) = i;   fj(end+1,1) = j;                      %#ok<AGROW>
        fx(end+1,1) = block_sum (Sa,i0,i1,j0,j1);                %#ok<AGROW>
        fy(end+1,1) = block_sum (Sb,i0,i1,j0,j1);                %#ok<AGROW>
      end
    end
  end
end
if isempty (fe), return; end
[sv, ord] = sort (fe, 'descend');
cs = cumsum (sv);
nt = find (cs >= prm.theta*cs(end), 1, 'first');
if isempty (nt), nt = numel (sv); end
mk = ord(1:nt);
big = max (fx, fy);
if prm.isotropic
  ex = true (size (big));  ey = ex;
else
  ex = fx >= prm.dir_kappa*big;
  ey = fy >= prm.dir_kappa*big;
  dead = ~ex & ~ey;  ex = ex | dead;  ey = ey | dead;
end
nm = {'U12','V1','V23','U3'};
newk = cell (L+1, 4);
for k = 1:numel (mk)
  m = mk(k);  l = fl(m);  q = fq(m);
  Pq = H.P(q);  brk = H.brk{l}{q};  nd = H.nd{l}{q};
  bu = brk{1}(:);  bv = brk{2}(:);
  i0 = max (1, fi(m)-p);  i1 = min (nd(1), fi(m));
  j0 = max (1, fj(m)-p);  j1 = min (nd(2), fj(m));
  if numel (H.Om) < l+1, H.Om{l+1} = zeros (0,4); end
  H.Om{l+1}(end+1,:) = [Pq.x0+Pq.Lx*bu(i0), Pq.x0+Pq.Lx*bu(i1+1), ...
                        Pq.y0+Pq.Ly*bv(j0), Pq.y0+Pq.Ly*bv(j1+1)];
  du = find (strcmp (nm, Pq.ku));   dv = find (strcmp (nm, Pq.kv));
  if ex(m)
    newk{l+1,du} = [newk{l+1,du}, 0.5*(bu(i0:i1).' + bu(i0+1:i1+1).')];
  end
  if ey(m)
    newk{l+1,dv} = [newk{l+1,dv}, 0.5*(bv(j0:j1).' + bv(j0+1:j1+1).')];
  end
end
Lnew = max (L, numel (H.Om));
for l = (L+1):Lnew
  H.KN{l} = H.KN{l-1};
  % do NOT reset H.Om{l} here: the loop above has just filled it with the
  % supports of the marked functions.
  if numel (H.Om) < l || isempty (H.Om{l}), H.Om{l} = zeros (0,4); end
end
for m = 2:numel (H.KN)                 % a knot must persist at deeper levels
  for d = 1:4
    add = [];
    for r = 2:min (m, size (newk,1))
      add = [add, newk{r,d}];          %#ok<AGROW>
    end
    H.KN{m}.(nm{d}) = insert_knots (H.KN{m}.(nm{d}), add, p);
  end
end
for l = 1:numel (H.KN)
  if numel (H.Om) < l, H.Om{l} = zeros (0,4); end
  H = build_level (H, l, prm);
end
H = rebuild_sets (H);
H = update_functions (H);
H = build_Csub (H);
end

function k = insert_knots (k, add, p)
if isempty (add), return; end
brk = unique (k);
add = unique (add(:).');
add = add(add > brk(1)+1e-13 & add < brk(end)-1e-13);
add = add(arrayfun (@(v) min (abs (v-brk)) > 1e-13, add));
if isempty (add), return; end
k = open_knots (sort ([brk, add]), p);
end

%%  TABLES

function print_tables (R, prm)
%PRINT_TABLES  iter/dof/levels/eta per (degree, driver); no exact solution, no I_eff.
ok = find ([R.ok]);
if isempty (ok), return; end
sep = repmat ('-', 1, 58);
tno = 0;
for k = ok
  tno = tno + 1;
  h = R(k).hist;  nd = h.ndof(:);  ev = h.eta(:);
  fprintf ('\n  %s\n', sep);
  fprintf ('  Table %d.  %s,  p = %d (C^%d),  eps = %.0e.\n', tno, ...
           driver_title (R(k).driver), R(k).degree, R(k).degree-1, prm.eps);
  fprintf ('  %s\n', sep);
  fprintf ('  %4s  %9s  %8s  %15s\n', 'k', 'dof', 'levels', eta_head (R(k).driver));
  fprintf ('  %s\n', sep);
  for i = 1:numel (nd)
    fprintf ('  %4d  %9d  %8d  %15.4e\n', i, nd(i), h.nlev(i), ev(i));
  end
  fprintf ('  %s\n', sep);
  fprintf (['  final eta = %.4e at %d dof (%d elements, %d levels);' ...
            '  rate vs sqrt(dof) = %.2f\n'], ev(end), nd(end), h.nel(end), ...
           h.nlev(end), -fit_slope_ (sqrt (nd), ev));
  fprintf ('  max aspect ratio %.2e\n', h.aspect(end));
end
fprintf ('\n');
end

function s = driver_title (drv)
if strcmpi (drv,'gradind'), s = 'gradient indicator (4.16)';
else,                       s = 'VMS residual estimator (4.13)';
end
end
function s = eta_head (drv)
if strcmpi (drv,'gradind'), s = 'eta_grad'; else, s = 'eta_vms'; end
end

function s = fit_slope_ (x, y)
lx = log (x(:));  ly = log (y(:));
v = isfinite (lx) & isfinite (ly) & (x(:) > 0) & (y(:) > 0);
if sum (v) < 2, s = NaN; return; end
cf = polyfit (lx(v), ly(v), 1);  s = cf(1);
end


%%  FIGURE INFRASTRUCTURE

function C = degree_colors ()
C = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.93 0.69 0.13];
end
function M = degree_markers ()
M = {'o','s','^','d','v'};
end
function ls = driver_style (drv)
if strcmpi (drv,'gradind'), ls = '--'; else, ls = '-'; end
end
function s = drv_name (drv)
if strcmpi (drv,'gradind'), s = '\eta_{grad}-driven';
else,                       s = '\eta_{vms}-driven';
end
end
function s = short_label (Rk)
if strcmpi (Rk.driver,'gradind'), s = sprintf ('p = %d,  \\eta_{grad}', Rk.degree);
else,                             s = sprintf ('p = %d,  \\eta_{vms}', Rk.degree);
end
end

function siam_axes (ax, prm)
set (ax, 'FontName', prm.fig_font, 'FontSize', prm.fig_fsize, ...
         'LineWidth', 0.75, 'TickDir', 'out', 'Box', 'on', ...
         'XGrid', 'on', 'YGrid', 'on', 'Layer', 'top');
try, set (ax, 'GridLineStyle', ':', 'GridAlpha', 0.30); catch, end
try, set (ax, 'MinorGridLineStyle','none','XMinorGrid','off','YMinorGrid','off'); catch, end
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
%SET_FIG_SIZE  Fixed pixel size, cascaded from previous figure, clamped on-screen.
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

function put_legend (ax, h, leg, key, ncol)
%PUT_LEGEND  Legend below axes, column-sorted by (degree, driver).
[~, ord] = sort (key);
lg = legend (ax, h(ord), leg(ord), 'Location', 'southoutside');
try, set (lg, 'NumColumns', max (1, ncol)); catch, end
try, set (lg, 'Box', 'off', 'FontSize', 10, 'Interpreter', 'tex'); catch, end
end

function save_figure (fig, name, prm)
if isempty (prm.fig_export), return; end
if ~exist (prm.fig_export, 'dir'), mkdir (prm.fig_export); end
try
  axs = findall (fig, 'Type', 'axes');
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
        case 'pdf', print (fig, '-dpdf', '-vector', fn);
        otherwise,  print (fig, '-dpng', '-r400',   fn);
      end
    end
    fprintf ('  figure written: %s\n', fn);
  catch werr
    warning ('thb:fig','could not write %s (%s).', fn, werr.message);
  end
end
end


%%  FIGURE 1: eta vs dof, every degree and both drivers on one axes

function plot_eta_all (R, prm)
ok = find ([R.ok]);
if isempty (ok), return; end
cascade_counter (true);
fig = figure ('Name','Estimator vs dof','Color','w');
set_fig_size (fig, 980, 660);
ax = axes ('Parent', fig);  hold (ax, 'on');
C = degree_colors ();  M = degree_markers ();
h = [];  leg = {};  key = [];  ndall = [];
for k = ok
  p = min (max (R(k).degree,1), 5);
  isg = strcmpi (R(k).driver, 'gradind');
  if isg, mfc = 'w'; else, mfc = C(p,:); end
  hh = plot (ax, R(k).hist.ndof, R(k).hist.eta, ...
             [M{p} driver_style(R(k).driver)], 'Color', C(p,:), ...
             'LineWidth', 1.6, 'MarkerSize', 6.0, ...
             'MarkerFaceColor', mfc, 'MarkerEdgeColor', C(p,:));
  h(end+1) = hh;  leg{end+1} = short_label (R(k));   %#ok<AGROW>
  key(end+1) = 10*R(k).degree + isg;                 %#ok<AGROW>
  ndall = [ndall; R(k).hist.ndof(:)];                %#ok<AGROW>
end
set (ax, 'XScale', 'log', 'YScale', 'log');
wide_log_xaxis (ax, ndall, prm.xpad);
xlabel (ax, 'degrees of freedom  N');
ylabel (ax, '\eta');
siam_axes (ax, prm);
put_legend (ax, h, leg, key, numel (unique ([R(ok).degree])));
hold (ax, 'off');
save_figure (fig, 'eta_vs_dof', prm);
end



%%  FIGURE SET 2: the final adaptive mesh, ONE FIGURE PER (degree, driver)

function plot_final_meshes (R, prm)
ok = find ([R.ok]);
if isempty (ok), return; end
cascade_counter (true);
for k = ok
  fig = figure ('Name', sprintf ('Final mesh, p=%d (%s)', R(k).degree, ...
                                 R(k).driver), 'Color','w');
  set_fig_size (fig, 560, 600);
  ax = axes ('Parent', fig);
  draw_mesh (ax, R(k).H);
  set (ax, 'FontName', prm.fig_font, 'FontSize', prm.fig_fsize);
  title (ax, sprintf ('p = %d, C^{%d},  %s,  %d elements,  %d dof', ...
                      R(k).degree, R(k).degree-1, drv_name (R(k).driver), ...
                      R(k).hist.nel(end), R(k).hist.ndof(end)));
  drawnow;
  save_figure (fig, sprintf ('final_mesh_p%d_%s', R(k).degree, ...
                             lower (R(k).driver)), prm);
end
end

function draw_mesh (ax, H)
%DRAW_MESH  Every ACTIVE cell of every level, so the hierarchy shows.
fill (ax, [0.5 1 1 0.5],[0 0 0.5 0.5],[0.93 0.93 0.93],'EdgeColor','none');
hold (ax,'on');
for l = 1:numel (H.KN)
  for q = 1:numel (H.P)
    if isempty (H.act{l}{q}), continue; end
    Pq = H.P(q);  brk = H.brk{l}{q};
    bx = Pq.x0 + Pq.Lx*brk{1}(:);   by = Pq.y0 + Pq.Ly*brk{2}(:);
    [i1, j1] = ind2sub (H.nd{l}{q}, H.act{l}{q}(:));
    X = [bx(i1).'; bx(i1+1).'; bx(i1+1).'; bx(i1).'];
    Y = [by(j1).'; by(j1).'; by(j1+1).'; by(j1+1).'];
    patch ('Parent',ax,'XData',X,'YData',Y,'FaceColor','none', ...
           'EdgeColor',[.25 .25 .25],'LineWidth',0.25);
  end
end
plot (ax,[0 0.5 0.5 1 1 0 0],[0 0 0.5 0.5 1 1 0],'k-','LineWidth',1.6);
plot (ax,0.5,0.5,'r^','MarkerSize',7,'MarkerFaceColor','r');
axis (ax,[0 1 0 1]); axis (ax,'square'); box (ax,'on');
set (ax,'XTick',[0 .5 1],'YTick',[0 .5 1],'Layer','top');
hold (ax,'off');
end


%%  FIGURE SET 3: the numerical solution, ONE FIGURE PER (degree, driver)

function plot_solution_surfs (R, prm)
ok = find ([R.ok]);
if isempty (ok), return; end
cascade_counter (true);
for k = ok
  H = R(k).H;  u = R(k).u;
  fig = figure ('Name', sprintf ('Solution, p=%d (%s)', R(k).degree, ...
                                 R(k).driver), 'Color','w');
  set_fig_size (fig, 640, 600);
  ax = axes ('Parent', fig);  hold (ax,'on');
  colormap (fig, parula);
  zmin = inf;  zmax = -inf;  drew = false;
  for q = 1:numel (H.P)
    try
      [X, Y, Z] = eval_patch_thb (H, u, q, prm.plot_npts);
    catch err
      warning ('thb:surf','patch %d: %s', q, err.message);  continue;
    end
    if isempty (Z) || ~any (isfinite (Z(:))), continue; end
    surf (ax, X, Y, Z, 'EdgeColor','none');
    zmin = min (zmin, min (Z(:)));  zmax = max (zmax, max (Z(:)));
    drew = true;
  end
  if ~drew, close (fig); continue; end
  shading (ax,'interp');
  xlabel (ax,'x'); ylabel (ax,'y'); zlabel (ax,'u_h');
  grid (ax,'on'); box (ax,'on'); view (ax, 42, 28);
  axis (ax,[0 1 0 1]);
  if zmax > zmin, zlim (ax, [zmin zmax]); end
  try, pbaspect (ax, [1 1 0.75]); catch, end
  set (ax,'FontName',prm.fig_font,'FontSize',prm.fig_fsize);
  title (ax, sprintf ('p = %d, C^{%d},  %s,  %d dof', R(k).degree, ...
                      R(k).degree-1, drv_name (R(k).driver), ...
                      R(k).hist.ndof(end)));
  hold (ax,'off');
  drawnow;
  save_figure (fig, sprintf ('solution_p%d_%s', R(k).degree, ...
                             lower (R(k).driver)), prm);
end
end

function [X, Y, Z] = eval_patch_thb (H, u, q, npts)
%EVAL_PATCH_THB  u_h on a tensor grid; each point evaluated at the finest active level.
p  = H.p;
ug = linspace (0,1,npts);  vg = linspace (0,1,npts);
Z  = nan (npts, npts);  todo = true (npts, npts);
for l = numel (H.KN) : -1 : 1
  if ~any (todo(:)), break; end
  if isempty (H.act{l}{q}), continue; end
  Nl  = size (H.C{l}, 2);
  uTP = H.C{l} * u(1:Nl);
  cp  = uTP(H.spmp{l}.gnum{q}(:));
  kn  = H.sp{l}{q}.knots;
  n1  = numel (kn{1}) - p - 1;
  n2  = numel (kn{2}) - p - 1;
  Cm  = reshape (cp, n1, n2);
  nd  = H.nd{l}{q};
  mask = false (nd);  mask(H.act{l}{q}) = true;
  brk = H.brk{l}{q};
  iu = cell_index_ (ug, brk{1}, nd(1));
  iv = cell_index_ (vg, brk{2}, nd(2));
  here = mask(iu, iv) & todo;
  if ~any (here(:)), continue; end
  [Bu, Cu] = bspline_values_ (kn{1}, p, ug);
  [Bv, Cv] = bspline_values_ (kn{2}, p, vg);
  Zl = zeros (npts, npts);
  for a = 1:npts
    T = Bu(a,:) * Cm(Cu(a,:), :);
    Zl(a,:) = sum (T(Cv) .* Bv, 2).';
  end
  Z(here) = Zl(here);
  todo(here) = false;
end
[U, V] = ndgrid (ug, vg);
X = H.P(q).x0 + H.P(q).Lx * U;
Y = H.P(q).y0 + H.P(q).Ly * V;
end

function i = cell_index_ (v, brk, n)
brk = brk(:).';
i = sum (bsxfun (@ge, v(:), brk(1:end-1)), 2);
i(v(:) >= brk(end)) = n;
i = min (max (i, 1), n);
end

function [B, cols] = bspline_values_ (knots, p, x)
knots = knots(:).';
n = numel (knots) - p - 1;
x = min (max (x(:).', knots(1)), knots(end));
s = findspan (n-1, p, x, knots);
B = basisfun (s, x, p, knots);
cols = bsxfun (@plus, double (s(:)) - p, 0:p) + 1;
cols = min (max (cols, 1), n);
end



%%  OPTIONS

function prm = parse_options (prm, args)
if numel (args) == 1 && isstruct (args{1})
  f = fieldnames (args{1});
  for i = 1:numel (f), prm.(f{i}) = args{1}.(f{i}); end
  return;
end
if mod (numel (args), 2) ~= 0
  error ('thb:args','Options must come in name/value pairs.');
end
for i = 1:2:numel (args)
  name = args{i};
  if ~ischar (name), error ('thb:args','Option names must be strings.'); end
  if isfield (prm, name),           prm.(name) = args{i+1};
  elseif isfield (prm.plots, name), prm.plots.(name) = args{i+1};
  else,  error ('thb:args','Unknown option ''%s''.', name);
  end
end
end