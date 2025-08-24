% ---------------------------------- %
% ------ORBITAL MECHANICS CODE------ %
% ---------------------------------- %

%{
Copyright (c) 2025 Sam Owen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
%}

clear; 
clc; 
close all;

%% ---------- Physical Constants and Nondimensional Scaling ----------
m1 = 5.974e24; % Earth mass [kg]
m2 = 7.348e22; % Moon mass [kg]
mu = m2 / (m1 + m2); % Nondimensional mass ratio (secondary)
R_em = 384400; % Earth-Moon distance [km]
Re = 6378.1363; % Earth radius [km]
moon_radius_km = 1737.4; % Moon radius [km]

%% ---------- User Parameters (Change These) ----------
N_particles = 1000; % Number of random test particles
altitude_km = 400; % Parking orbit altitude above earth surface [km]
max_inclination_deg = 0; % Maximum random inclination

dv_frac_max = 0.45; % Max |deltaV| as fraction of local circular speed
dv_frac_min = 0.35; % Min |deltaV| as fraction of local circular speed
sim_time_nd = 6; % Nondimensional integration time
npoints = 3000; % Output points per integration
capture_thresh_km = 5 * moon_radius_km; % Success threshold

burn = char('xy'); % Set 'r' for fully random burn, 'xyz' for prograde weighted with z-component, 'xy' for prograde weighted without z-component

%% ---------- Precompute Useful Quantities ----------
r_orbit_km = Re + altitude_km; % Orbit radius from Earth's centre [km]
r_orbit_n  = r_orbit_km / R_em; % Nondimensional

% Approximate local circular speed (nondimensional) about Earth (two-body approx)
v_circ_n = sqrt((1 - mu) / r_orbit_n);

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
tt_eval = linspace(0, sim_time_nd, npoints);

% Storage
trajectories = cell(N_particles,1);
results = struct('min_dist_km', cell(N_particles,1), 'reached', cell(N_particles,1), 'dv_frac', cell(N_particles,1), 'inc_deg', cell(N_particles,1), 'raan_deg', cell(N_particles,1), 'nu_deg', cell(N_particles,1), 'dv_dir', cell(N_particles, 1));

rng('shuffle'); % Randomise

moon_pos_rot_km = [(1-mu)*R_em, 0, 0]; % Moon in rotating frame (dimensional)
G = 6.67430e-20; % Gravitational constant [km^3/kg/s^2]
m_tot = m1 + m2;
omega = sqrt(G * m_tot / R_em^3); % [rad/s]
T = 1 / omega; % Characteristic time scale [s] (4.34 days)

tic;

%% ---------- Monte-Carlo Particle Generation & Single-Impulse ----------
for k = 1:N_particles
    % Random true anomaly and orbital orientation
    nu = 2*pi*rand(); % True anomaly [rad]
    inc = deg2rad(max_inclination_deg * rand()); % Inclination [rad] 
    raan = 2*pi * rand(); % RAAN [rad]

    % Position in perifocal frame (nondimensional)
    r_pf = r_orbit_n * [cos(nu); sin(nu); 0];
    v_pf = v_circ_n * [-sin(nu); cos(nu); 0]; % Circular velocity in perifocal frame

    % Rotation: R = Rz(raan) * Rx(inc) (arg of perigee set to zero)
    Rz_raan = [cos(raan) -sin(raan) 0; sin(raan) cos(raan) 0; 0 0 1];
    Rx_inc  = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
    R_orb = Rz_raan * Rx_inc;

    % Inertial (Earth-centred) position & velocity (nondimensional)
    r_eci = R_orb * r_pf;
    v_eci = R_orb * v_pf;

    % Convert to barycentric rotating-frame coordinates: Earth centre is at x=-mu
    r_bary = [-mu; 0; 0] + r_eci; % Nondim

    % Rotating-frame velocity: v_rot = v_inertial - Omega x r = v_eci + [y; -x; 0]
    v_rot = v_eci + [r_bary(2); -r_bary(1); 0];

    % Get local velocity unit vector in rotating frame
    vel_dir = v_rot / norm(v_rot);

    if strcmpi(burn, 'r')
        dv_dir = randn(3,1); dv_dir = dv_dir / norm(dv_dir);

    elseif strcmpi(burn, 'xyz')
        % Pick a weighted mix of prograde and random
        bias_strength = 0.8; % 0 = random, 1 = fully prograde 
        dv_dir = bias_strength * vel_dir + (1 - bias_strength) * randn(3,1);

    elseif strcmpi(burn, 'xy')
        % Pick a weighted mix of prograde and random
        bias_strength = 0.8; % 0 = random, 1 = fully prograde 
        rand_xy = randn(2,1); 
        dv_dir_xy = bias_strength * vel_dir(1:2) + (1 - bias_strength) * rand_xy; 
        dv_dir_xy = dv_dir_xy / norm(dv_dir_xy); 
        dv_dir = [dv_dir_xy; 0]; 
    end

    dv_frac = dv_frac_min + rand() * (dv_frac_max - dv_frac_min);
    deltaV = dv_dir * (dv_frac * v_circ_n);

    v_rot_post = v_rot + deltaV;

    % Initial state for CR3BP integration
    Y0 = [r_bary; v_rot_post];

    % Integrate
    [tt, Y] = ode45(@(t,Y) nondim_cr3bp(t, Y, mu), tt_eval, Y0, opts);

    % Compute dimensional rotating-frame positions [km] and distance to Moon
    r_rot_km = Y(:,1:3) * R_em;
    d2moon = sqrt( sum( (r_rot_km - moon_pos_rot_km).^2 , 2) );
    [min_d_km, min_idx] = min(d2moon);
    reached = (min_d_km <= capture_thresh_km);

    % Store
    trajectories{k}.tt = tt;
    trajectories{k}.Y = Y;
    trajectories{k}.min_d_km = min_d_km;
    trajectories{k}.min_idx = min_idx;
    trajectories{k}.reached = reached;
    trajectories{k}.dv_frac = dv_frac;
    trajectories{k}.inc_deg = rad2deg(inc);
    trajectories{k}.raan_deg = rad2deg(raan);
    trajectories{k}.nu_deg = rad2deg(nu);

    results(k).min_dist_km = min_d_km;
    results(k).reached = reached;
    results(k).dv_frac = dv_frac;
    results(k).inc_deg = rad2deg(inc);
    results(k).raan_deg = rad2deg(raan);
    results(k).nu_deg = rad2deg(nu);
    results(k).dv_dir = dv_dir;

    clc;
    fprintf('Calculating...\n');
    fprintf('Iteration Number %.0f / %.0f', k, N_particles);
end

clc;

calc_time = toc;
fprintf('Simulation Time: %.2f seconds\n', calc_time);
fprintf('Lunar Capture Radius: %.2f km\n', capture_thresh_km);

%% ---------- Report & Plot Results (Rotating Frame) ----------
success_idx = find([results.reached]);
fprintf('\nRan %d Trials - Found %d Successes\n\n', N_particles, numel(success_idx));

if ~isempty(success_idx)
    for i = 1:numel(success_idx)
        k = success_idx(i);
        fprintf('#%d: True Anomaly = %.1f deg, Inclination = %.1f deg, RAAN = %.1f deg, DeltaV = %.3f km/s\n', ...
            i, results(k).nu_deg, results(k).inc_deg, results(k).raan_deg, results(k).dv_frac * v_circ_n * (R_em / T));
    end
else
    fprintf('No Close Approaches Found\n');
end

%% ---------- Time ----------
hill_radius_km = R_em * (m2/(3*m1))^(1/3); % Lunar SOI radius [km]
fprintf('\n');

for i = 1:numel(success_idx)
    k = success_idx(i);

    if results(k).reached
        r_rot_km = trajectories{k}.Y(:,1:3) * R_em; 
        r_rel = r_rot_km - moon_pos_rot_km; 
        r_mag = sqrt(sum(r_rel.^2, 2)); 
        t_dim = tt * T;  

        % Find first time inside capture threshold
        idx_capture = find(r_mag <= capture_thresh_km, 1, 'first');

        if ~isempty(idx_capture)
            time_to_capture = t_dim(idx_capture);

            % Find all indices inside the SOI
            inside_soi_idx = find(r_mag <= hill_radius_km);

            if ~isempty(inside_soi_idx)
                % Break into continuous segments
                inside_diff = diff(inside_soi_idx);
                break_points = find(inside_diff > 1);

                segment_starts = [inside_soi_idx(1); inside_soi_idx(break_points+1)];
                segment_ends   = [inside_soi_idx(break_points); inside_soi_idx(end)];

                % Sum durations
                total_time = sum(t_dim(segment_ends) - t_dim(segment_starts));
            else
                total_time = 0;
            end

            fprintf('#%d: Time to Lunar Capture Radius = %.2f hours, Total Time in Lunar SOI = %.2f hours\n', ...
                i, time_to_capture/3600, total_time/3600);
        end
    end
end

%% ---------- First Figure: All Trajectories ----------
figure('Color','k','Position',[200 200 900 700]);
ax = axes; hold on; axis equal; grid on; box on;
ax.Color = 'k'; ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';
xlabel('x [km]','Color','w'); ylabel('y [km]','Color','w'); zlabel('z [km]','Color','w');

view(3);

% Plot primaries
% Create unit sphere
[sx, sy, sz] = sphere(50); 

% Scale spheres by actual radii
earth_radius = Re; % Earth radius [km]
moon_radius = moon_radius_km; % Moon radius [km]

% Earth position in rotating frame (x = -mu * R_em)
earth_center = [-mu * R_em, 0, 0];
% Moon position in rotating frame (x = (1 - mu) * R_em)
moon_center = [(1 - mu) * R_em, 0, 0];

hold on;

% Plot Earth
surf(earth_center(1) + earth_radius * sx, ...
     earth_center(2) + earth_radius * sy, ...
     earth_center(3) + earth_radius * sz, ...
     'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'Earth');

% Plot Moon
surf(moon_center(1) + moon_radius * sx, ...
     moon_center(2) + moon_radius * sy, ...
     moon_center(3) + moon_radius * sz, ...
     'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', 'Moon'); 

% Plot trajectories
for k = 1:N_particles
    rrot = trajectories{k}.Y(:,1:3) * R_em;
    if trajectories{k}.reached
        plot3(rrot(:,1), rrot(:,2), rrot(:,3), 'g', 'LineWidth', 1.2);
    else
        plot3(rrot(:,1), rrot(:,2), rrot(:,3), 'r', 'LineWidth', 0.6);
    end
end

%% ---------- Second Figure: Only Successful Trajectories ----------
figure('Color','k','Position',[300 200 900 700]);
ax2 = axes; hold on; axis equal; grid on; box on;
ax2.Color = 'k'; ax2.XColor = 'w'; ax2.YColor = 'w'; ax2.ZColor = 'w';
xlabel('x [km]','Color','w'); ylabel('y [km]','Color','w'); zlabel('z [km]','Color','w');

view(3);

% Plot primaries
% Plot Earth
surf(earth_center(1) + earth_radius * sx, ...
     earth_center(2) + earth_radius * sy, ...
     earth_center(3) + earth_radius * sz, ...
     'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'Earth');

% Plot Moon
surf(moon_center(1) + moon_radius * sx, ...
     moon_center(2) + moon_radius * sy, ...
     moon_center(3) + moon_radius * sz, ...
     'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', 'Moon'); 

% Define some colors for successful trajectories
colors = lines(numel(success_idx)); 

for i = 1:numel(success_idx)
    k = success_idx(i);
    rrot = trajectories{k}.Y(:,1:3) * R_em;
    plot3(rrot(:,1), rrot(:,2), rrot(:,3), 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Trajectory #%d', i));
end

legend('TextColor', 'w');

%% ---------- Third Figure: Animated Trajectories ----------
figure('Color','k','Position',[300 200 900 700]);
ax2 = axes; hold on; axis equal; grid on; box on;
ax2.Color = 'k'; ax2.XColor = 'w'; ax2.YColor = 'w'; ax2.ZColor = 'w';
xlabel('x [km]','Color','w'); ylabel('y [km]','Color','w'); zlabel('z [km]','Color','w');

view(3);

% Plot primaries
% Plot Earth
surf(earth_center(1) + earth_radius * sx, ...
     earth_center(2) + earth_radius * sy, ...
     earth_center(3) + earth_radius * sz, ...
     'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'Earth');

% Plot Moon
surf(moon_center(1) + moon_radius * sx, ...
     moon_center(2) + moon_radius * sy, ...
     moon_center(3) + moon_radius * sz, ...
     'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', 'Moon'); 

colors = lines(numel(success_idx)); 

moon_x = moon_center(1);

animatedLines = gobjects(numel(success_idx),1);
for i = 1:numel(success_idx)
    animatedLines(i) = plot3(nan, nan, nan, 'Color', colors(i,:), 'LineWidth', 1.5);
end

moon_threshold_km = capture_thresh_km;

max_len = max(cellfun(@(t) size(t.Y,1), trajectories(success_idx)));

finished = false(numel(success_idx),1);

for idx = 1:max_len
    for i = 1:numel(success_idx)
        if finished(i)
            continue; % Skip finished trajectories
        end
        k = success_idx(i);
        rrot = trajectories{k}.Y(:,1:3) * R_em;

        if idx <= size(rrot,1)

            set(animatedLines(i), 'XData', rrot(1:idx,1), 'YData', rrot(1:idx,2), 'ZData', rrot(1:idx,3));
            
            dist_to_moon = sqrt( (rrot(idx,1) - moon_x)^2 + rrot(idx,2)^2 + rrot(idx,3)^2 );
            if dist_to_moon <= moon_threshold_km
                finished(i) = true;
            end
        else
            finished(i) = true; 
        end
    end
    
    drawnow;
    pause(0.02); % Adjust animation speed here
    
    if all(finished)
        break;
    end
end

%% ---------- CR3BP EoM ----------
function Ydot = nondim_cr3bp(~, Y, mu)
    x = Y(1); y = Y(2); z = Y(3);
    xdot = Y(4); ydot = Y(5); zdot = Y(6);
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    Ydot = zeros(6,1);
    Ydot(1:3) = [xdot; ydot; zdot];
    Ydot(4) = 2*ydot + x - (1 - mu)*(x + mu)/r1^3 - mu*(x - 1 + mu)/r2^3;
    Ydot(5) = -2*xdot + y - (1 - mu)*y/r1^3 - mu*y/r2^3;
    Ydot(6) = -(1 - mu)*z/r1^3 - mu*z/r2^3;
end