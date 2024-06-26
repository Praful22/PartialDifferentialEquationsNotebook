% Parameters
L = 1;          % Length of string
T = 30;         % Total simulation time
Nx = 100;       % Number of spatial points
Nt = 10000;     % Number of time steps (increased for stability)
c = 1;          % Wave speed
zeta = 0.05;    % Damping ratio
omega_0 = pi;   % Natural frequency
k_air = 0.01;   % Air resistance coefficient (reduced)

dx = L / (Nx-1);
dt = T / (Nt-1);
x = linspace(0, L, Nx);
t = linspace(0, T, Nt);

% Check CFL condition
CFL = c * dt / dx;
if CFL > 1
    error('CFL condition not met. Decrease dt or increase dx.');
end

% Initial condition
f = @(x) 0.1 * sin(pi*x/L);  % Initial shape of the string

% Initialize solution array
u = zeros(Nx, Nt);
u(:,1) = f(x);
u(:,2) = u(:,1);  % Initial velocity is zero

% Time-stepping (using a more stable scheme)
for n = 2:Nt-1
    % Interior points
    u(2:Nx-1,n+1) = (2-2*CFL^2-2*zeta*omega_0*dt)*u(2:Nx-1,n) + ...
                    (CFL^2-k_air*dt*abs(u(2:Nx-1,n)-u(2:Nx-1,n-1))/dt).*(u(3:Nx,n)+u(1:Nx-2,n)) + ...
                    (-1+2*zeta*omega_0*dt+k_air*dt*abs(u(2:Nx-1,n)-u(2:Nx-1,n-1))/dt).*u(2:Nx-1,n-1) - ...
                    omega_0^2*dt^2*u(2:Nx-1,n);
    
    % Boundary conditions
    u(1,n+1) = 0;  % Fixed at x = 0
    u(Nx,n+1) = 0; % Fixed at x = L
end

% Create animation
fig = figure;
for n = 1:100:Nt  % Reduced number of frames for faster animation
    plot(x, u(:,n), 'LineWidth', 2);
    axis([0 L -0.15 0.15]);
    title(sprintf('Time: %.2f', t(n)));
    xlabel('Position');
    ylabel('Displacement');
    drawnow;
    
    % Capture the plot as an image
    frame = getframe(fig);
    im{n} = frame2im(frame);
end

% Save animation as GIF
filename = 'string_vibration_realistic_damping.gif';
for n = 1:100:Nt
    [A,map] = rgb2ind(im{n},256);
    if n == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end

% Create and save surface plot
figure;
[T, X] = meshgrid(t, x);
surf(T, X, u);
shading interp;
colormap(jet);
colorbar;
xlabel('Time');
ylabel('Position');
zlabel('Displacement');
title('String Vibration with Realistic Damping');
saveas(gcf, 'string_vibration_surface_realistic_damping.png');

% Plot maximum displacement over time
figure;
plot(t, max(abs(u), [], 1));
xlabel('Time');
ylabel('Maximum Displacement');
title('Damping Effect on String Motion');
saveas(gcf, 'max_displacement_over_time_realistic_damping.png');