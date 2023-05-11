clc, clear, close all

% parameters
mu = 398600.5;
J2 = 0.00108263;
R = 6378 ;           % Equatorial radius (R/r <1)

% Initial Conditions from Textbook's Example 10.6
R0 = [-2384.46; 5729.01; 3050.46];   % [km]
V0 = [-7.36138; -2.98997; 1.64354];  % [km/s]
X0 = [R0; V0];
[a, e, i0, omega0, w0, f0, h0] = rv2coe(R0, V0, mu);

% Confirm that the following CoEs are obtained
a = 8059;
h0 = 55839;
e = 0.17136;
f0 = 40 *pi/180;
omega0 = 45 *pi/180;
i0 = 28 * pi/180;
w0 = 30 * pi/180;

C0 = [ h0, e,  f0, omega0, i0,  w0 ];

% Determine Orbital Period
n = sqrt(mu/a^3);
T = 2*pi/n;

%% Set-up

% Simulation Time
sim_time1 = 48*60*60;    % second
sim_time2 = 5*T;
sim_time3 = 30*60;   % 1 hour

% Erorr tolerance option for ODE45
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Unperturbed orbit ODE
EoM = @(t, x)[zeros(3,3), eye(3);...
    -(mu/norm(x(:,1))^3)*eye(3), zeros(3,3)]*x;
[time, State] = ode45(EoM, [0 sim_time1], X0, options);

% J2 pertubed orbit ODE from Gauss Variational Equation
[time_p, State_p] = ode45(@gauss_var, [0 sim_time1], C0, options);

% pre-allocate memory
r = zeros(3, length(time_p));                 
v = zeros(3, length(time_p)); 
rx = zeros(1, length(time_p));
ry = zeros(1, length(time_p));
rz = zeros(1, length(time_p));
vx = zeros(1, length(time_p));
vy = zeros(1, length(time_p));
vz = zeros(1, length(time_p));

a_p = zeros(1, length(time_p));  % semi-major axis

%----------------------------------Animation ----------------------------------------------
% curve = animatedline('LineWidth', 1);
% set(gca, 'XLim', [-1e4 1e4], 'YLim', [-1e4 1e4], 'ZLim', [-1e4, 4e4]);
% grid on, hold on
% view(-42,42)
% title('Molniya Orbit','FontSize',16)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % head = scatter3(Earth(1), Earth(2), Earth(3), 2000, 'filled', 'MarkerFaceColor','b')
% set(gcf, 'WindowState', 'maximized');
%---------------------------------------------------------------------------------------

%% Computing position and velocity from CoE
for j = 1:length(time_p)

    h = State_p(j, 1);  omega = State_p(j, 4);
    e = State_p(j, 2);      i = State_p(j, 5);
    f = State_p(j, 3);      w = State_p(j, 6);
    a_p(j) = (h^2)/(mu*(1-e^2));      % semi-major axis
    
    r = h^2/(mu*(1 + e*cos(f)));
    
    r_rsw = [r; 0; 0];
    v_rsw = [h * e * sin(f)/(a*(1 - e^2)); h/r; 0];
    
    [r_ijk, v_ijk] = rsw2ijk(r_rsw, v_rsw, omega, i, w, f);
    
    rx(j) = r_ijk(1);
    ry(j) = r_ijk(2);
    rz(j) = r_ijk(3);
    
    vx(j) = v_ijk(1);
    vy(j) = v_ijk(2);
    vz(j) = v_ijk(3);
    
    %---------------------------Animation-------------------------------------------------
%     title(num2str(time_p(j)/(60*60),'time = %4.0f hrs'));
%     addpoints(curve, r_ijk(1), r_ijk(2), r_ijk(3));
%     head = scatter3(r_ijk(1), r_ijk(2), r_ijk(3), 50, 'filled', 'MarkerFaceColor','r');
%     drawnow
%     delete(head)
    %-------------------------------------------------------------------------------------
end

%% Plot CoE Variations
LineWidth = 1.5;

figure(1)
subplot(6,1,1)
plot(time_p/(60*60), (State_p(:,4)-State_p(1,4))*180/pi, 'LineWidth',LineWidth)
title('RAAN Variations (deg)')

subplot(6,1,2)
plot(time_p/(60*60), (State_p(:,6)-State_p(1,6))*180/pi, 'LineWidth',LineWidth)
title('Argument of Perigee Variations (deg)')

subplot(6,1,3)
plot(time_p/(60*60), (State_p(:,1)-State_p(1,1)), 'LineWidth',LineWidth)
title('Angular Momentum Variations (km^2/s)')

subplot(6,1,4)
plot(time_p/(60*60), State_p(:,2)-State_p(1,2), 'LineWidth',LineWidth)
title('Eccentricity Variations')

subplot(6,1,5)
plot(time_p/(60*60), (State_p(:,5)-State_p(1,5)) *180/pi, 'LineWidth',LineWidth)
title('Inclination Variations (deg)')

subplot(6,1,6)
plot(time_p/(60*60), a_p-a, 'LineWidth',LineWidth)
title('Semi-major Axis (km)')
xlabel('Time (hours)')

%% Plot the Orbit
figure('units','normalized','outerposition',[0 0 1 1])
% plot perturbed orbit
plot3(rx, ry, rz, '-b', 'LineWidth', 1.5)
hold on, axis equal

% plot unperturbed orbit
plot3(State(:,1), State(:,2), State(:,3),'LineWidth', 1.5)

% plot Earth
[x,y,z] = sphere;
x = x*R;   % Earth real radius
y = y*R;
z = z*R;

surf(x,y,z)
% shading interp      % turn this on to remove Earth's edges
alpha 0.25
colormap parula
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off

title('J2 perturbed orbit','FontSize',16)
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on
hold off

legend("J2 Perturbed Orbit", "Unperturbed Orbit", "Earth")