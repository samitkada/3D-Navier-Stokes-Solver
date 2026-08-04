%% ============================================================
%  3D LID-DRIVEN CUBIC CAVITY FLOW
%  Solves incompressible Navier-Stokes equations in 3D
%  using finite differences + pressure-velocity coupling (SIMPLE-like)
%  
%  Physics: A cube of fluid. The top lid moves at U=1 in the x-direction.
%           Viscosity resists the motion. A 3D vortex structure forms.
%
%  Tools: MATLAB | Author: [Samit Kadasinghanahalli]
%% ============================================================

clc; clear; close all;

%% ============================================================
%  SECTION 1: PARAMETERS
%  Start with a coarse grid (N=16) to verify it runs, then increase
%% ============================================================

N    = 32;          % grid points per side — try 16 first, then 24, then 32
Re   = 400;         % Reynolds number — try 100 (smooth) up to 400 (complex)
                    % Re = U*L/nu,  here U=1, L=1, so nu = 1/Re
nu   = 1 / Re;      % kinematic viscosity

U_lid = 1.0;        % lid velocity in x-direction (top face, z=1)

dt    = 0.0005;      % time step — reduce if it diverges
n_steps = 8000;     % total time steps — increase for higher Re
save_every = 100;   % how often to print progress

%% ============================================================
%  SECTION 2: GRID SETUP
%  Uniform Cartesian grid on unit cube [0,1]^3
%  Staggered grid: u,v,w at cell faces, p at cell centers
%  For simplicity here we use a collocated grid
%% ============================================================

h  = 1 / (N - 1);              % grid spacing
x  = linspace(0, 1, N);
y  = linspace(0, 1, N);
z  = linspace(0, 1, N);

[X, Y, Z] = meshgrid(x, y, z);  % [N x N x N] grids

% Initialize velocity and pressure fields
u = zeros(N, N, N);   % x-velocity
v = zeros(N, N, N);   % y-velocity
w = zeros(N, N, N);   % z-velocity
p = zeros(N, N, N);   % pressure

% Apply lid BC: top face (k=N) moves at U_lid in x
u(:, :, N) = U_lid;

fprintf('\n========================================\n');
fprintf('  3D LID-DRIVEN CAVITY SOLVER\n');
fprintf('========================================\n');
fprintf('Grid:      %d x %d x %d\n', N, N, N);
fprintf('Re:        %d\n', Re);
fprintf('nu:        %.4f\n', nu);
fprintf('dt:        %.4f\n', dt);
fprintf('Steps:     %d\n', n_steps);
fprintf('========================================\n');
fprintf('Running... (this may take a few minutes)\n\n');

%% ============================================================
%  SECTION 3: TIME INTEGRATION
%  Fractional step method (Chorin projection):
%    1. Compute intermediate velocity u* (ignore pressure)
%    2. Solve Poisson equation for pressure correction
%    3. Project u* onto divergence-free field
%
%  Spatial derivatives: 2nd-order central differences
%% ============================================================

tic;
for step = 1:n_steps

    %% -- Step 3a: Compute intermediate velocities (u*, v*, w*) --
    % Interior indices only (boundary conditions fixed)
    i = 2:N-1;  j = 2:N-1;  k = 2:N-1;

    % Convective terms (upwind for stability)
    % u-momentum
    dudx = (u(i,j+1,k) - u(i,j-1,k)) / (2*h);
    dudy = (u(i+1,j,k) - u(i-1,j,k)) / (2*h);
    dudz = (u(i,j,k+1) - u(i,j,k-1)) / (2*h);

    % v-momentum
    dvdx = (v(i,j+1,k) - v(i,j-1,k)) / (2*h);
    dvdy = (v(i+1,j,k) - v(i-1,j,k)) / (2*h);
    dvdz = (v(i,j,k+1) - v(i,j,k-1)) / (2*h);

    % w-momentum
    dwdx = (w(i,j+1,k) - w(i,j-1,k)) / (2*h);
    dwdy = (w(i+1,j,k) - w(i-1,j,k)) / (2*h);
    dwdz = (w(i,j,k+1) - w(i,j,k-1)) / (2*h);

    % Diffusion terms (Laplacian)
    lap_u = (u(i,j+1,k) + u(i,j-1,k) + u(i+1,j,k) + u(i-1,j,k) ...
           + u(i,j,k+1) + u(i,j,k-1) - 6*u(i,j,k)) / h^2;

    lap_v = (v(i,j+1,k) + v(i,j-1,k) + v(i+1,j,k) + v(i-1,j,k) ...
           + v(i,j,k+1) + v(i,j,k-1) - 6*v(i,j,k)) / h^2;

    lap_w = (w(i,j+1,k) + w(i,j-1,k) + w(i+1,j,k) + w(i-1,j,k) ...
           + w(i,j,k+1) + w(i,j,k-1) - 6*w(i,j,k)) / h^2;

    % Pressure gradients
    dpdx = (p(i,j+1,k) - p(i,j-1,k)) / (2*h);
    dpdy = (p(i+1,j,k) - p(i-1,j,k)) / (2*h);
    dpdz = (p(i,j,k+1) - p(i,j,k-1)) / (2*h);

    % Intermediate velocity (explicit Euler, no pressure yet)
    u_star = u;
    v_star = v;
    w_star = w;

    u_star(i,j,k) = u(i,j,k) + dt * (-u(i,j,k).*dudx ...
                                      -v(i,j,k).*dudy ...
                                      -w(i,j,k).*dudz ...
                                      + nu*lap_u - dpdx);

    v_star(i,j,k) = v(i,j,k) + dt * (-u(i,j,k).*dvdx ...
                                      -v(i,j,k).*dvdy ...
                                      -w(i,j,k).*dvdz ...
                                      + nu*lap_v - dpdy);

    w_star(i,j,k) = w(i,j,k) + dt * (-u(i,j,k).*dwdx ...
                                      -v(i,j,k).*dwdy ...
                                      -w(i,j,k).*dwdz ...
                                      + nu*lap_w - dpdz);

    % Re-apply boundary conditions on u*
    u_star(:,:,N) = U_lid;   % lid
    u_star(:,:,1) = 0;       % bottom
    u_star(:,1,:) = 0;       % front
    u_star(:,N,:) = 0;       % back
    u_star(1,:,:) = 0;       % left
    u_star(N,:,:) = 0;       % right
    v_star([1 N],:,:) = 0;
    v_star(:,[1 N],:) = 0;
    v_star(:,:,[1 N]) = 0;
    w_star([1 N],:,:) = 0;
    w_star(:,[1 N],:) = 0;
    w_star(:,:,[1 N]) = 0;

    %% -- Step 3b: Solve pressure Poisson equation --
    % div(u*) = RHS for pressure correction
    div_u = (u_star(i,j+1,k) - u_star(i,j-1,k)) / (2*h) ...
           +(v_star(i+1,j,k) - v_star(i-1,j,k)) / (2*h) ...
           +(w_star(i,j,k+1) - w_star(i,j,k-1)) / (2*h);

    % Iterative pressure solve (Gauss-Seidel, 20 iterations)
    p_new = p;
    for piter = 1:20
        p_new(i,j,k) = (p_new(i,j+1,k) + p_new(i,j-1,k) ...
                      + p_new(i+1,j,k) + p_new(i-1,j,k) ...
                      + p_new(i,j,k+1) + p_new(i,j,k-1) ...
                      - h^2 * div_u / dt) / 6;
        % Neumann BC on pressure (dp/dn = 0 at walls)
        p_new(1,:,:)   = p_new(2,:,:);
        p_new(N,:,:)   = p_new(N-1,:,:);
        p_new(:,1,:)   = p_new(:,2,:);
        p_new(:,N,:)   = p_new(:,N-1,:);
        p_new(:,:,1)   = p_new(:,:,2);
        p_new(:,:,N)   = p_new(:,:,N-1);
    end
    p = p_new;

    %% -- Step 3c: Project velocity (correct with pressure) --
    u(i,j,k) = u_star(i,j,k) - dt*(p(i,j+1,k) - p(i,j-1,k))/(2*h);
    v(i,j,k) = v_star(i,j,k) - dt*(p(i+1,j,k) - p(i-1,j,k))/(2*h);
    w(i,j,k) = w_star(i,j,k) - dt*(p(i,j,k+1) - p(i,j,k-1))/(2*h);

    % Re-apply BCs after projection
    u(:,:,N) = U_lid;
    u(:,:,1) = 0; u(:,1,:) = 0; u(:,N,:) = 0;
    u(1,:,:) = 0; u(N,:,:) = 0;
    v(:)     = v .* (1 - (X==0|X==1|Y==0|Y==1|Z==0|Z==1));
    w(:)     = w .* (1 - (X==0|X==1|Y==0|Y==1|Z==0|Z==1));

    %% -- Progress report --
    if mod(step, save_every) == 0
        vel_mag = sqrt(u.^2 + v.^2 + w.^2);
        fprintf('Step %4d / %d  |  max|V| = %.4f  |  t = %.3f s\n', ...
                step, n_steps, max(vel_mag(:)), step*dt);
    end
end

elapsed = toc;
fprintf('\nDone! Elapsed: %.1f seconds\n\n', elapsed);

%% ============================================================
%  SECTION 4: POST-PROCESSING & 3D VISUALIZATION
%% ============================================================

vel_mag = sqrt(u.^2 + v.^2 + w.^2);

% Vorticity components
[dudy_f, dudx_f, dudz_f] = gradient(u, h);
[dvdy_f, dvdx_f, dvdz_f] = gradient(v, h);
[dwdy_f, dwdx_f, dwdz_f] = gradient(w, h);

omega_x = dwdy_f - dvdz_f;   % x-vorticity
omega_y = dudz_f - dwdx_f;   % y-vorticity
omega_z = dvdx_f - dudy_f;   % z-vorticity
omega_mag = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);

%% ============================================================
%  SECTION 5: PLOTS
%% ============================================================

fig = figure('Name', '3D Lid-Driven Cavity', ...
             'Color', [0.06 0.06 0.10], ...
             'Position', [60 40 1300 860]);

ax_bg = [0.10 0.10 0.16];
style_ax = @(ax) set(ax, 'Color', ax_bg, ...
    'XColor',[0.7 0.7 0.7],'YColor',[0.7 0.7 0.7],'ZColor',[0.7 0.7 0.7],...
    'GridColor',[0.22 0.22 0.30],'GridAlpha',1,...
    'XGrid','on','YGrid','on','ZGrid','on',...
    'FontSize',9,'FontName','Consolas');

%% --- Plot 1: Mid-plane velocity magnitude (z=0.5 slice) ---
ax1 = subplot(2,3,1);
k_mid = round(N/2);
contourf(squeeze(X(:,:,k_mid)), squeeze(Y(:,:,k_mid)), ...
         squeeze(vel_mag(:,:,k_mid)), 20, 'LineColor','none');
hold on;
% Velocity vectors (subsample)
step_q = max(1, round(N/10));
quiver(squeeze(X(1:step_q:end,1:step_q:end,k_mid)), ...
       squeeze(Y(1:step_q:end,1:step_q:end,k_mid)), ...
       squeeze(u(1:step_q:end,1:step_q:end,k_mid)), ...
       squeeze(v(1:step_q:end,1:step_q:end,k_mid)), ...
       0.5, 'w', 'LineWidth', 0.8);
colormap(ax1, jet); colorbar;
style_ax(ax1); axis equal tight;
xlabel('x'); ylabel('y');
title('Velocity  |  z=0.5 slice', 'Color','w','FontSize',11);

%% --- Plot 2: Vorticity magnitude mid-plane ---
ax2 = subplot(2,3,2);
contourf(squeeze(X(:,:,k_mid)), squeeze(Y(:,:,k_mid)), ...
         squeeze(omega_mag(:,:,k_mid)), 20, 'LineColor','none');
colormap(ax2, hot); colorbar;
style_ax(ax2); axis equal tight;
xlabel('x'); ylabel('y');
title('Vorticity magnitude  |  z=0.5', 'Color','w','FontSize',11);

%% --- Plot 3: u-velocity profile along vertical centerline ---
ax3 = subplot(2,3,3);
j_mid = round(N/2);
plot(squeeze(u(:,j_mid,k_mid)), y, ...
     'Color',[0.27 0.82 1.0],'LineWidth',2); hold on;
xline(0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.8);
style_ax(ax3);
xlabel('u-velocity'); ylabel('y');
title('U profile at centerline', 'Color','w','FontSize',11);
text(0.02, 0.08, sprintf('Re = %d', Re), ...
     'Color','w','FontSize',10,'FontName','Consolas','Units','normalized');

%% --- Plot 4: 3D velocity magnitude isosurface ---
ax4 = subplot(2,3,[4 5 6]);
iso_val = 0.3 * max(vel_mag(:));
p_iso = patch(isosurface(X, Y, Z, vel_mag, iso_val));
p_iso.FaceColor = [0.27 0.82 1.0];
p_iso.EdgeColor = 'none';
p_iso.FaceAlpha = 0.35;
hold on;

% Second isosurface (higher speed core)
iso_val2 = 0.6 * max(vel_mag(:));
p_iso2 = patch(isosurface(X, Y, Z, vel_mag, iso_val2));
p_iso2.FaceColor = [1.0 0.45 0.15];
p_iso2.EdgeColor = 'none';
p_iso2.FaceAlpha = 0.6;

% Vorticity isosurface
iso_om = 0.4 * max(omega_mag(:));
p_om = patch(isosurface(X, Y, Z, omega_mag, iso_om));
p_om.FaceColor = [1.0 0.85 0.2];
p_om.EdgeColor = 'none';
p_om.FaceAlpha = 0.2;

% Lighting for 3D depth
camlight('headlight');
lighting gouraud;

% Cavity outline
plot3([0 1 1 0 0],[0 0 0 0 0],[0 0 1 1 0],'w-','LineWidth',0.8);
plot3([0 1 1 0 0],[1 1 1 1 1],[0 0 1 1 0],'w-','LineWidth',0.8);
plot3([0 0],[0 1],[0 0],'w-','LineWidth',0.8);
plot3([1 1],[0 1],[0 0],'w-','LineWidth',0.8);
plot3([0 0],[0 1],[1 1],'w-','LineWidth',0.8);
plot3([1 1],[0 1],[1 1],'w-','LineWidth',0.8);

% Lid indicator
patch([0 1 1 0],[0 0 1 1],[1 1 1 1], ...
      'FaceColor',[1 1 1],'FaceAlpha',0.08,'EdgeColor','w','LineWidth',1);
text(0.5, 0.5, 1.05, sprintf('Lid  U=%.0f →', U_lid), ...
     'Color','w','FontSize',9,'HorizontalAlignment','center','FontName','Consolas');

style_ax(ax4);
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('3D Isosurfaces  |  Velocity (blue/orange) & Vorticity (yellow)  |  Re=%d', Re), ...
      'Color','w','FontSize',11);
view(-38, 22);
legend([p_iso p_iso2 p_om], ...
       {'Low speed', 'High speed', 'Vortex core'}, ...
       'TextColor','w','Color',ax_bg,'EdgeColor',[0.3 0.3 0.4],'FontSize',8);

sgtitle(sprintf('3D Lid-Driven Cavity  |  N=%d  |  Re=%d  |  t=%.2fs', ...
        N, Re, n_steps*dt), ...
        'Color','w','FontSize',14,'FontWeight','bold','FontName','Consolas');