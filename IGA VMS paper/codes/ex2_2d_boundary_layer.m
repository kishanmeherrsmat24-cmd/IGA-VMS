%% IGA SGS/SUPG -- VMS residual estimator, uniform then adaptive refinement
%
%      -eps*Lap(u) + b.grad(u) + c*u = f in (0,1)^2,  u = u_ex on dOmega
%      b = (2,3)^T,  c = 1
%   Runs UNIFORM then ADAPTIVE (Doerfler on VMS eta_K) for each p in p_list.
%
%   tau_stab = (h/(2|b|))(coth alpha - 1/alpha),  alpha = |b|h/(2 eps)
%   tau_up   = min( CA*h/|b| , CD*h^2/eps ),  CA=1/sqrt(2), CD=1/(3*sqrt(10))
%   eta_K    = tau_up,K * ||R||_{L2(K)},  I_eff = eta / ||u-u_h||_{L2}
%
%   Requires GeoPDEs 3.x and the NURBS toolbox.

%% USER PARAMETERS 
b_val   = [2; 3];          % convection field
c_val   = 1;               % reaction coefficient
b_norm  = norm(b_val);

p_list   = [1,2,3,4,5];    % B-spline degrees (C^{p-1})
eps_list = [1e-4];         % diffusion coefficients
nel_list = [4,8,16,32,64]; % uniform meshes, elements per direction
eps_fig  = 1e-4;           % which eps the summary figures show

% --- adaptive settings ---------------------------------------------------
CFG.nel0     = 4;          % initial adaptive mesh
CFG.theta    = 0.40;       % Dorfler bulk fraction
CFG.n_adapt  = 40;         % max adaptive iterations
CFG.dir_kappa= 0.50;       % bisection-direction test, see mark_directions
CFG.grade_ratio = 2;       % neighbouring spans within this factor (Inf = off)

% --- stopping controls ---------------------------------------------------
CFG.tol_rel  = 1e-13;
CFG.tol_abs  = 0;
CFG.max_ndof = 5000;
CFG.h_min    = 0;
CFG.stall_max= 3;

% --- stabilisation -------------------------------------------------------
CFG.h_stab   = 'min_edge';    % 'max_chord'|'flow_extent'|'stream_harm'|
                              % 'mean_chord'|'geometric'|'max_edge'|'min_edge'
CFG.test_op  = 'sgs';         % 'sgs' | 'supg' | 'gls'
CFG.stabilise= true;          % false -> plain Galerkin
CFG.p_corr   = false;         % true -> divide tau_stab by p

% --- VMS estimator -------------------------------------------------------
CFG.CA       = 1/sqrt(2);     % advective constant
CFG.CD       = 1/(3*sqrt(10));% diffusive constant
CFG.h_adv    = 'max_chord';   % advective branch of tau_up
CFG.h_dif    = 'max_chord';   % diffusive branch of tau_up
CFG.tau_rule = 'min';         % 'min' | 'harmonic'

% --- quadrature / housekeeping -------------------------------------------
CFG.nquad_add = 3;         % Gauss points per direction = p + nquad_add
CFG.chunk_mem = 1e6;       % elements evaluated in chunks of this size
CFG.b        = b_val;
CFG.c        = c_val;
CFG.bnorm    = b_norm;
CFG.ahat     = b_val / b_norm;


%% Exact solution and RHS
make_pde = @(ev) struct( ...
'eps',   ev, ...
'u_ex',  @(x,y) x.*y.^2 - y.^2.*exp(2*(x-1)/ev) ...
                  - x.*exp(3*(y-1)/ev) + exp((2*(x-1)+3*(y-1))/ev), ...
'f_rhs', @(x,y) x.*y.^2 + 6*x.*y + 2*y.^2 - 2*ev*x ...
                  + (2*ev - y.^2 - 6*y).*exp(2*(x-1)/ev) ...
                  - (x + 2).*exp(3*(y-1)/ev) ...
                  + exp((2*(x-1)+3*(y-1))/ev) );

if ~any(abs(eps_list - eps_fig) < 1e-14*max(1,abs(eps_list)))
    error('iga:epsfig','eps_fig = %g is not one of eps_list.', eps_fig);
end

%% Storage for summary figures (one entry per p, at eps = eps_fig)
np = numel(p_list);
sum_uni = cell(np,1);
sum_ada = cell(np,1);

%%  loop over p
for pi = 1:np
    p = p_list(pi);
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('   p = %d   (B-spline degree, C^%d)\n', p, p-1);
    fprintf('================================================================\n');

    %% --- (A) Uniform refinement table -----------------------------------
    for ei = 1:numel(eps_list)
        ev  = eps_list(ei);
        pde = make_pde(ev);
        nk  = numel(nel_list);
        ndof = zeros(1,nk);  etav = zeros(1,nk);  errs = zeros(1,nk);  iev = zeros(1,nk);
        for ki = 1:nk
            nrb            = build_uniform_nrb(nel_list(ki), p);
            [sp, msh, brk] = build_mesh_space(nrb, p, p + CFG.nquad_add);
            u_h            = solve_stab(sp, msh, brk, p, pde, CFG);
            E              = vms_estimate(sp, msh, brk, u_h, p, pde, CFG);
            ndof(ki) = sp.ndof;
            etav(ki) = E.eta;
            errs(ki) = E.err;
            iev(ki)  = E.eta / max(E.err, realmin);
        end
        fprintf('\n  Table.  Uniform refinement,  p = %d,  eps = %.0e\n', p, ev);
        print_iter_table(1:nk, ndof, etav, errs, iev);
        if abs(ev - eps_fig) < 1e-14*max(1,ev)
            U.ndof = ndof;  U.l2 = errs;  U.eta = etav;  U.ieff = iev;
            sum_uni{pi} = U;
        end
    end

    %% --- (B) Adaptive refinement, VMS-driven ----------------------------
    fprintf('\n  Adaptive SGS/SUPG, VMS-driven  (p=%d, theta=%.2f, h_stab=''%s'', test_op=''%s'')\n', ...
            p, CFG.theta, CFG.h_stab, CFG.test_op);
    for ei = 1:numel(eps_list)
        ev  = eps_list(ei);
        pde = make_pde(ev);
        fprintf('\n  eps = %.0e\n', ev);
        H = run_adaptive(p, pde, CFG);
        if abs(ev - eps_fig) < 1e-14*max(1,ev)
            sum_ada{pi} = H;
        end
    end
    fprintf('\n');
end

%%  Summary figures -- L2 error and effectivity index, uniform vs adaptive
Cc = degree_colors();   Mk = degree_markers();

figure('Position',[60 80 900 680],'Color','w');
ax1 = axes; hold(ax1,'on');
for pi = 1:np
    p  = p_list(pi);
    loglog(ax1, sum_uni{pi}.ndof(:), sum_uni{pi}.l2(:), [Mk{p} '-'], ...
        'Color', Cc(p,:), 'LineWidth',1.8,'MarkerSize',6.5,'MarkerFaceColor',Cc(p,:), ...
        'DisplayName', sprintf('p=%d, C^{%d}', p, p-1));
end
set(ax1,'XScale','log','YScale','log','FontSize',13,'LineWidth',1.2);
grid(ax1,'on'); grid(ax1,'minor'); box(ax1,'on');
xlabel(ax1,'degrees of freedom  N','FontSize',15,'FontWeight','bold');
ylabel(ax1,'L^2 error','FontSize',15,'FontWeight','bold');
legend(ax1,'Location','southwest','FontSize',9,'Box','on');
hold(ax1,'off');

figure('Position',[100 80 900 680],'Color','w');
ax2 = axes; hold(ax2,'on');
for pi = 1:np
    p = p_list(pi);
    H = sum_ada{pi};
    loglog(ax2, H.ndof(:), H.l2(:), [Mk{p} '-'], 'Color', Cc(p,:), ...
        'LineWidth',1.9,'MarkerSize',7,'MarkerFaceColor',Cc(p,:), ...
        'DisplayName', sprintf('p=%d, C^{%d}', p, p-1));
end
set(ax2,'XScale','log','YScale','log','FontSize',13,'LineWidth',1.2);
grid(ax2,'on'); grid(ax2,'minor'); box(ax2,'on');
xlabel(ax2,'degrees of freedom  N','FontSize',15,'FontWeight','bold');
ylabel(ax2,'L^2 error','FontSize',15,'FontWeight','bold');
legend(ax2,'Location','southwest','FontSize',9,'Box','on');
hold(ax2,'off');

figure('Position',[140 60 900 680],'Color','w');
ax3 = axes; hold(ax3,'on');
for pi = 1:np
    p = p_list(pi);
    loglog(ax3, sum_uni{pi}.ndof(:), sum_uni{pi}.ieff(:), [Mk{p} '-'], ...
        'Color', Cc(p,:), 'LineWidth',1.8,'MarkerSize',6.5,'MarkerFaceColor',Cc(p,:), ...
        'DisplayName', sprintf('p=%d, C^{%d}', p, p-1));
end
set(ax3,'XScale','log','YScale','log');
xl = xlim(ax3);  yl = ylim(ax3);
patch(ax3, [xl(1) xl(2) xl(2) xl(1)], [1 1 10 10], [0.85 0.92 0.85], ...
'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
plot(ax3, xl, [1 1], 'k-', 'LineWidth',1.2,'HandleVisibility','off');
try, uistack(findobj(ax3,'Type','patch'),'bottom'); catch, end
xlim(ax3, xl);  ylim(ax3, yl);
grid(ax3,'on'); grid(ax3,'minor'); box(ax3,'on');
set(ax3,'FontSize',13,'LineWidth',1.2);
xlabel(ax3,'degrees of freedom  N','FontSize',15,'FontWeight','bold');
ylabel(ax3,'I_{eff} = \eta / || u - u_h ||_{L^2}','FontSize',15,'FontWeight','bold');
legend(ax3,'Location','eastoutside','FontSize',9,'Box','on');
hold(ax3,'off');

figure('Position',[180 40 900 680],'Color','w');
ax4 = axes; hold(ax4,'on');
for pi = 1:np
    p = p_list(pi);
    H = sum_ada{pi};
    loglog(ax4, H.ndof(:), H.ieff(:), [Mk{p} '-'], 'Color', Cc(p,:), ...
        'LineWidth',1.9,'MarkerSize',7,'MarkerFaceColor',Cc(p,:), ...
        'DisplayName', sprintf('p=%d, C^{%d}', p, p-1));
end
set(ax4,'XScale','log','YScale','log');
xl = xlim(ax4);  yl = ylim(ax4);
patch(ax4, [xl(1) xl(2) xl(2) xl(1)], [1 1 10 10], [0.85 0.92 0.85], ...
'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
plot(ax4, xl, [1 1], 'k-', 'LineWidth',1.2,'HandleVisibility','off');
try, uistack(findobj(ax4,'Type','patch'),'bottom'); catch, end
xlim(ax4, xl);  ylim(ax4, yl);
grid(ax4,'on'); grid(ax4,'minor'); box(ax4,'on');
set(ax4,'FontSize',13,'LineWidth',1.2);
xlabel(ax4,'degrees of freedom  N','FontSize',15,'FontWeight','bold');
ylabel(ax4,'I_{eff} = \eta / || u - u_h ||_{L^2}','FontSize',15,'FontWeight','bold');
legend(ax4,'Location','eastoutside','FontSize',9,'Box','on');
hold(ax4,'off');

%%  FINAL SUMMARY  (degree, run, dof, eta, L2, Ieff -- final iteration only)
final_summary(sum_uni, sum_ada, p_list, eps_fig);

%%  LOCAL FUNCTIONS

% ---------------------------------------------------------------- geometry
function nrb = build_uniform_nrb(nel, p)
    kv0   = [0 0 1 1];
    coefs = zeros(2,2,2);
    coefs(1,2,1)=1;  coefs(1,2,2)=1;
    coefs(2,1,2)=1;  coefs(2,2,2)=1;
    ik  = linspace(0,1,nel+1);  ik = ik(2:end-1);
    nrb = nrbmak(coefs,{kv0,kv0});
    nrb = nrbdegelev(nrb,[p-1,p-1]);
if ~isempty(ik), nrb = nrbkntins(nrb,{ik,ik}); end
end

function [sp, msh, brk] = build_mesh_space(nrb, p, nquad)
%BUILD_MESH_SPACE  No sp_precompute: shape functions evaluated in chunks
%   (full hessian ~1GB at p=5, N=64).
    geo     = geo_load(nrb);
    brk     = {unique(nrb.knots{1}), unique(nrb.knots{2})};
    [qn,qw] = msh_set_quad_nodes(brk, msh_gauss_nodes([nquad, nquad]));
    msh     = msh_cartesian(brk, qn, qw, geo);
    sp      = sp_bspline(nrb.knots, [p,p], msh);
end

% ----------------------------------------------------------- element sizes
function [hx, hy, ix, iy] = elem_geom(brk, msh_el, iel)
%ELEM_GEOM  True knot-span widths and (i,j) span index of element IEL.
    cx = mean(msh_el.geo_map(1,:,iel));
    cy = mean(msh_el.geo_map(2,:,iel));
    ix = find(brk{1} <= cx + 1e-13, 1, 'last');
    iy = find(brk{2} <= cy + 1e-13, 1, 'last');
    ix = min(max(ix,1), numel(brk{1})-1);
    iy = min(max(iy,1), numel(brk{2})-1);
    hx = brk{1}(ix+1) - brk{1}(ix);
    hy = brk{2}(iy+1) - brk{2}(iy);
end

function h = elem_length(hx, hy, ahat, mode)
%ELEM_LENGTH  Scalar element length for the chosen mode.
    a1 = abs(ahat(1));  a2 = abs(ahat(2));  tiny = 1e-14;
switch lower(mode)
case 'max_chord',   h = min(hx / max(a1,tiny), hy / max(a2,tiny));
case 'mean_chord',  h = hx * hy / max(a2*hx + a1*hy, tiny);
case 'flow_extent', h = a1*hx + a2*hy;
case 'stream_harm', h = 2 / (a1/hx + a2/hy);
case 'geometric',   h = sqrt(hx * hy);
case 'max_edge',    h = max(hx, hy);
case 'min_edge',    h = min(hx, hy);
otherwise, error('iga:hmode','Unknown element-length mode ''%s''.', mode);
end
end

function v = xi_coth(a)
% coth(a) - 1/a, with correct a->0 limit a/3
if a < 1e-8, v = a/3; else, v = coth(a) - 1/a; end
end

function tau = tau_stab(hx, hy, p, pde, CFG)
% tau_stab = (h/(2 pd |b|))(coth alpha - 1/alpha),  alpha = |b|h/(2 pd eps)
    hf  = elem_length(hx, hy, CFG.ahat, CFG.h_stab);
    pd  = 1;  if CFG.p_corr, pd = p; end
    al  = CFG.bnorm * hf / (2 * pd * pde.eps);
    tau = hf / (2 * pd * CFG.bnorm) * xi_coth(al);
end

function t = tau_up(hx, hy, pde, CFG)
% tau_up = min( CA*h_adv/|b| , CD*h_dif^2/eps )
    ha  = elem_length(hx, hy, CFG.ahat, CFG.h_adv);
    hd  = elem_length(hx, hy, CFG.ahat, CFG.h_dif);
    adv = max(CFG.CA * ha / CFG.bnorm, realmin);
    dif = max(CFG.CD * hd^2 / pde.eps, realmin);
switch lower(CFG.tau_rule)
case 'min',      t = min(adv, dif);
case 'harmonic', t = 1 / (1/adv + 1/dif);
otherwise, error('iga:taurule','Unknown tau_rule ''%s''.', CFG.tau_rule);
end
end

% ------------------------------------------------------------ chunked eval
function lst = elem_chunks(msh, sp, p, CFG)
    nsh = (p+1)^2;
    nc  = max(1, floor(CFG.chunk_mem / max(msh.nqn * nsh, 1)));
    nc  = min(nc, msh.nel);
    lst = cell(1, ceil(msh.nel/nc));
for k = 1:numel(lst)
        lst{k} = ((k-1)*nc + 1) : min(k*nc, msh.nel);
end
end

function [msh_el, sp_el] = eval_chunk(msh, sp, elems, want_hess)
if nargin < 4, want_hess = false; end
    msh_el = msh_evaluate_element_list(msh, elems);
if want_hess
        sp_el = sp_evaluate_element_list(sp, msh_el, 'value', true, ...
'gradient', true, 'hessian', true);
        H = sp_el.shape_function_hessians;
        sp_el.shape_function_laplacians = ...
             reshape(H(1,1,:,:,:) + H(2,2,:,:,:), size(H,3), size(H,4), size(H,5));
return;
end
try
        sp_el = sp_evaluate_element_list(sp, msh_el, 'value', true, ...
'gradient', true, 'laplacian', true);
catch
        sp_el = sp_evaluate_element_list(sp, msh_el, 'value', true, ...
'gradient', true, 'hessian', true);
        H = sp_el.shape_function_hessians;
        sp_el.shape_function_laplacians = ...
             reshape(H(1,1,:,:,:) + H(2,2,:,:,:), size(H,3), size(H,4), size(H,5));
end
end

function [N, Gx, Gy, Lp] = shape_blocks(sp_el, iel, nqn, nsh)
    N  = reshape(sp_el.shape_functions(:,:,iel),              nqn, nsh);
    Gx = reshape(sp_el.shape_function_gradients(1,:,:,iel),   nqn, nsh);
    Gy = reshape(sp_el.shape_function_gradients(2,:,:,iel),   nqn, nsh);
    Lp = reshape(sp_el.shape_function_laplacians(:,:,iel),    nqn, nsh);
end

% ----------------------------------------------------------------- solver
function u = solve_stab(sp, msh, brk, p, pde, CFG)
%SOLVE_STAB  Galerkin + stabilisation, chunked triplet assembly.
    b = CFG.b;  cc = CFG.c;  ev = pde.eps;
    nsh = (p+1)^2;
    I = zeros(msh.nel*nsh*nsh, 1);  J = I;  V = I;  ptr = 0;
    F = zeros(sp.ndof,1);
for ch = elem_chunks(msh, sp, p, CFG)
        elems = ch{1};
        [msh_el, sp_el] = eval_chunk(msh, sp, elems);
        nqn = msh_el.nqn;  nshm = sp_el.nsh_max;
for iel = 1:msh_el.nel
            w  = msh_el.quad_weights(:,iel) .* msh_el.jacdet(:,iel);
            x  = reshape(msh_el.geo_map(1,:,iel), nqn, 1);
            y  = reshape(msh_el.geo_map(2,:,iel), nqn, 1);
            [N, Gx, Gy, Lp] = shape_blocks(sp_el, iel, nqn, nshm);
            bG = b(1)*Gx + b(2)*Gy;
            R  = -ev*Lp + bG + cc*N;
if CFG.stabilise
                [hx, hy] = elem_geom(brk, msh_el, iel);
                tau = tau_stab(hx, hy, p, pde, CFG);
switch lower(CFG.test_op)
case 'sgs',  P = bG + ev*Lp - cc*N;
case 'supg', P = bG;
case 'gls',  P = R;
otherwise, error('iga:testop','Unknown test_op ''%s''.', CFG.test_op);
end
else
                tau = 0;  P = zeros(nqn, nshm);
end
            fq = pde.f_rhs(x, y);
            W  = w(:, ones(1, nshm));
            Ke = N.'*(W.*bG) + ev*(Gx.'*(W.*Gx) + Gy.'*(W.*Gy)) ...
                 + cc*(N.'*(W.*N)) + tau*(P.'*(W.*R));
            Fe = N.'*(w.*fq) + tau*(P.'*(w.*fq));
            conn = sp_el.connectivity(:,iel);
            idx  = ptr + (1:nshm*nshm);
            I(idx) = repmat(conn, nshm, 1);
            J(idx) = reshape(repmat(conn.', nshm, 1), [], 1);
            V(idx) = Ke(:);
            ptr = ptr + nshm*nshm;
            F(conn) = F(conn) + Fe;
end
end
    A = sparse(I(1:ptr), J(1:ptr), V(1:ptr), sp.ndof, sp.ndof);
    % Dirichlet data by L2 projection of u_ex on each boundary edge
    bc_dofs_all = [];  u_bc_all = [];
for k = 1:4
        bmsh  = msh_precompute(msh.boundary(k));
        bsp   = sp_precompute(sp.boundary(k), bmsh, 'value',true,'gradient',false);
        M_bnd = op_u_v(bsp, bsp, bmsh, ones(bmsh.nqn, bmsh.nel));
        F_bnd = op_f_v(bsp, bmsh, pde.u_ex(reshape(bmsh.geo_map(1,:,:), bmsh.nqn, bmsh.nel), ...
                                           reshape(bmsh.geo_map(2,:,:), bmsh.nqn, bmsh.nel)));
        bc_dofs_all = [bc_dofs_all; bsp.dofs(:)];   %#ok<AGROW>
        u_bc_all    = [u_bc_all;    M_bnd\F_bnd];   %#ok<AGROW>
end
    [bc_dofs, ia] = unique(bc_dofs_all, 'last');
    u_bc     = u_bc_all(ia);
    int_dofs = setdiff(1:sp.ndof, bc_dofs);
    u = zeros(sp.ndof,1);
    u(bc_dofs)  = u_bc;
    u(int_dofs) = A(int_dofs,int_dofs) \ ...
                  (F(int_dofs) - A(int_dofs,bc_dofs)*u_bc);
end

% -------------------------------------------------------------- estimator
function E = vms_estimate(sp, msh, brk, u, p, pde, CFG)
%VMS_ESTIMATE  Per-element eta_K, true L2 error, and Hessian direction weights.
%   E.eta_K = tau_up,K * ||R_h||_{L2(K)}
    b = CFG.b;  cc = CFG.c;  ev = pde.eps;
    n = msh.nel;
    eta_K  = zeros(n,1);  err_K  = zeros(n,1);
    ax_K   = zeros(n,1);  ay_K   = zeros(n,1);
    ix_K   = zeros(n,1);  iy_K   = zeros(n,1);
for ch = elem_chunks(msh, sp, p, CFG)
        elems = ch{1};
        [msh_el, sp_el] = eval_chunk(msh, sp, elems, true);
        H = sp_el.shape_function_hessians;
        nqn = msh_el.nqn;  nshm = sp_el.nsh_max;
for iel = 1:msh_el.nel
            k  = elems(iel);
            w  = msh_el.quad_weights(:,iel) .* msh_el.jacdet(:,iel);
            x  = reshape(msh_el.geo_map(1,:,iel), nqn, 1);
            y  = reshape(msh_el.geo_map(2,:,iel), nqn, 1);
            [N, Gx, Gy, Lp] = shape_blocks(sp_el, iel, nqn, nshm);
            ue = u(sp_el.connectivity(:,iel));
            uh = N*ue;  gx = Gx*ue;  gy = Gy*ue;  lp = Lp*ue;
            Rh = -ev*lp + b(1)*gx + b(2)*gy + cc*uh - pde.f_rhs(x, y);
            [hx, hy, ix, iy] = elem_geom(brk, msh_el, iel);
            eta_K(k)  = tau_up(hx, hy, pde, CFG) * sqrt(max(sum(w.*Rh.^2), 0));
            err_K(k) = sqrt(max(sum(w.*(pde.u_ex(x,y) - uh).^2), 0));
            ix_K(k)  = ix;  iy_K(k) = iy;
            % Hessian-based directional weights for anisotropic marking
            Hxx = reshape(H(1,1,:,:,iel), nqn, nshm);
            Hyy = reshape(H(2,2,:,:,iel), nqn, nshm);
            ax_K(k) = hx^2 * sqrt(max(sum(w.*(Hxx*ue).^2), 0));
            ay_K(k) = hy^2 * sqrt(max(sum(w.*(Hyy*ue).^2), 0));
end
end
    E.eta_K = eta_K;  E.err_K = err_K;
    E.ax = ax_K;      E.ay = ay_K;
    E.ix = ix_K;      E.iy = iy_K;
    E.eta  = sqrt(sum(eta_K.^2));
    E.err  = sqrt(sum(err_K.^2));
end

% ------------------------------------------------------- adaptive driver
function H = run_adaptive(p, pde, CFG)
%RUN_ADAPTIVE  SOLVE -> ESTIMATE -> MARK -> REFINE loop.
    nrb = build_uniform_nrb(CFG.nel0, p);

    n = CFG.n_adapt;
    H.ndof = zeros(n,1);  H.l2 = zeros(n,1);  H.eta = zeros(n,1);  H.ieff = zeros(n,1);

    it_done = 0;  ind_ref = NaN;  ind_best = Inf;  stall = 0;
    why = 'iteration limit';

    fmt_sep = repmat('-', 1, 46);
    fprintf('  %s\n', fmt_sep);
    fprintf('  %4s  %8s  %12s  %12s  %8s\n', 'it', 'dof', 'eta', 'L2', 'Ieff');
    fprintf('  %s\n', fmt_sep);

for iter = 1:n
        [sp, msh, brk] = build_mesh_space(nrb, p, p + CFG.nquad_add);
        u_h = solve_stab(sp, msh, brk, p, pde, CFG);
        E   = vms_estimate(sp, msh, brk, u_h, p, pde, CFG);

        H.ndof(iter) = sp.ndof;
        H.l2(iter)   = E.err;
        H.eta(iter)  = E.eta;
        H.ieff(iter) = E.eta / max(E.err, realmin);
        it_done      = iter;

        fprintf('  %4d  %8d  %12.4e  %12.4e  %8.3f\n', ...
                iter, sp.ndof, E.eta, E.err, H.ieff(iter));

        ind_K = E.eta_K;
        ind   = E.eta;
if iter == 1, ind_ref = max(ind, 1e-300); end

        hx = diff(brk{1}(:));  hy = diff(brk{2}(:));

        % Stopping controls
if ind >= ind_best, stall = stall + 1; else, stall = 0; end
        ind_best = min(ind_best, ind);
if ind / ind_ref < CFG.tol_rel, why = 'relative indicator tolerance'; break; end
if ind < CFG.tol_abs,           why = 'absolute indicator tolerance'; break; end
if stall >= CFG.stall_max,      why = 'indicator stopped decreasing'; break; end
if sp.ndof > CFG.max_ndof,      why = 'dof budget';                  break; end
if iter == n, break; end

        % Mark, grade, refine
        mk = dorfler_mark(ind_K, CFG.theta);
        [mx, my] = mark_directions(mk, E, CFG.dir_kappa, numel(hx), numel(hy));
        mx = mx & (hx(:).' > 2*CFG.h_min);
        my = my & (hy(:).' > 2*CFG.h_min);
        mx = grade_marks(brk{1}, mx, CFG.grade_ratio);
        my = grade_marks(brk{2}, my, CFG.grade_ratio);
        mx = mx & (hx(:).' > 2*CFG.h_min);
        my = my & (hy(:).' > 2*CFG.h_min);

        nrb_new = insert_bisect(nrb, brk, mx, my);
if isequal(nrb_new.knots{1}, nrb.knots{1}) && ...
           isequal(nrb_new.knots{2}, nrb.knots{2})
            why = 'no span left to bisect (h_min reached)';  break;
end
        nrb = nrb_new;
end

    f = {'ndof','l2','eta','ieff'};
for q = 1:numel(f), H.(f{q}) = H.(f{q})(1:it_done); end
    H.it_done = it_done;  H.why = why;

    fprintf('  %s\n', fmt_sep);
    fprintf('  Stopped after %d iterations: %s\n', it_done, why);
    fprintf('  %s\n', fmt_sep);
end

% ------------------------------------------------------ marking / refining
function marked = dorfler_mark(eta_K, theta)
%DORFLER_MARK  Bulk marking on the VMS indicator.
    e2 = eta_K(:).^2;
if all(e2 == 0), marked = false(size(e2)); return; end
    [~, idx] = sort(e2, 'descend');
    cs  = cumsum(e2(idx));
    n_m = find(cs >= theta*sum(e2), 1, 'first');
if isempty(n_m), n_m = numel(e2); end
    marked = false(numel(e2),1);
    marked(idx(1:max(n_m,1))) = true;
end

function [mx, my] = mark_directions(marked, E, kappa, nspx, nspy)
%MARK_DIRECTIONS  Which direction to bisect based on Hessian content.
%   Split in x if hx^2||d2u/dx2||_K >= kappa*max(ax,ay), likewise in y.
    m  = marked(:);
    ax = E.ax(:);  ay = E.ay(:);
    big = max(ax, ay);
    ex  = m & (ax >= kappa * big);
    ey  = m & (ay >= kappa * big);
    dead = m & ~ex & ~ey;    % both weights zero -> split both ways
    ex = ex | dead;   ey = ey | dead;
    mx = false(1, nspx);  my = false(1, nspy);
    mx(E.ix(ex)) = true;
    my(E.iy(ey)) = true;
end

function m = grade_marks(brk, m, ratio)
%GRADE_MARKS  Extend marks until bisected spans stay within RATIO of neighbours.
    m = logical(m(:)).';
if ~isfinite(ratio) || ratio <= 1 || ~any(m), return; end
    h0 = diff(brk(:)).';
for sweep = 1:100
        h   = h0 ./ (1 + double(m));
        add = false(size(m));
        add(1:end-1) = add(1:end-1) | (h(1:end-1) > ratio*h(2:end)   & ~m(1:end-1));
        add(2:end)   = add(2:end)   | (h(2:end)   > ratio*h(1:end-1) & ~m(2:end));
if ~any(add), return; end
        m = m | add;
end
end

function nrb_new = insert_bisect(nrb, brk, mx, my)
%INSERT_BISECT  New knot at midpoint of every marked span.
    tol = 1e-13;
    new_kx = midpoints(brk{1}, mx, unique(nrb.knots{1}), tol);
    new_ky = midpoints(brk{2}, my, unique(nrb.knots{2}), tol);
    nrb_new = nrb;
if ~isempty(new_kx), nrb_new = nrbkntins(nrb_new, {new_kx, []}); end
if ~isempty(new_ky), nrb_new = nrbkntins(nrb_new, {[], new_ky}); end
end

function nk = midpoints(brk_d, mask, existing, tol)
    nk = [];
    idx = find(logical(mask(:)).');
if isempty(idx), return; end
    idx = idx(idx >= 1 & idx <= numel(brk_d)-1);
if isempty(idx), return; end
    nk = 0.5*(brk_d(idx) + brk_d(idx+1));
    nk = unique(nk(:).');
    nk = nk(nk > tol & nk < 1-tol);
if ~isempty(nk)
        nk = nk(arrayfun(@(v) min(abs(v - existing)) > tol, nk));
end
end

% -------------------------------------------------------- table / summary
function print_iter_table(it, ndof, eta, l2, ieff)
sep = repmat('-', 1, 46);
fprintf('  %s\n', sep);
fprintf('  %4s  %8s  %12s  %12s  %8s\n', 'it', 'dof', 'eta', 'L2', 'Ieff');
fprintf('  %s\n', sep);
for k = 1:numel(it)
    fprintf('  %4d  %8d  %12.4e  %12.4e  %8.3f\n', it(k), ndof(k), eta(k), l2(k), ieff(k));
end
fprintf('  %s\n', sep);
end

function final_summary(sum_uni, sum_ada, p_list, eps_fig)
%FINAL_SUMMARY  One table row per (p, run): degree, dof, eta, L2, Ieff.
    np = numel(p_list);
    wid = repmat('=', 1, 74);
    fprintf('\n\n  %s\n', wid);
    fprintf('   FINAL SUMMARY -- VMS residual estimator      (eps = %.0e)\n', eps_fig);
    fprintf('  %s\n', wid);
    fprintf('   %-2s | %-9s | %8s | %12s | %12s | %8s\n', ...
'p', 'run', 'dof', 'eta', 'L2', 'Ieff');
    fprintf('  %s\n', repmat('-', 1, 74));
for pi = 1:np
        p = p_list(pi);
        U = sum_uni{pi};  A = sum_ada{pi};
        fprintf('   %2d | %-9s | %8d | %12.4e | %12.4e | %8.3f\n', ...
                p, 'uniform', U.ndof(end), U.eta(end), U.l2(end), U.ieff(end));
        fprintf('   %2s | %-9s | %8d | %12.4e | %12.4e | %8.3f\n', ...
'', 'adaptive', A.ndof(end), A.eta(end), A.l2(end), A.ieff(end));
if pi < np, fprintf('  %s\n', repmat('-', 1, 74)); end
end
    fprintf('  %s\n\n', wid);
end

% ------------------------------------------------------------------ plots
function C = degree_colors()
C = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; ...
     0.49 0.18 0.56; 0.93 0.69 0.13];
end

function M = degree_markers()
M = {'o','s','^','d','v'};
end