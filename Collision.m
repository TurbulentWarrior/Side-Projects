clear;
clc;

containerSize=5;
N_obj = 4;
obj_size = 0.1*ones(N_obj);  % The size represents the radius of the objects
obj_mass = ones(N_obj);

% Create coordinates
obj_x = zeros(N_obj, 1);
obj_y = zeros(N_obj, 1);

% Create velocities
obj_vx = zeros(N_obj, 1);
obj_vy = zeros(N_obj, 1);

% Initialize coordinates and velocities
obj_x(4) = 0.75;
obj_y(4) = 2;
obj_x(3) = 1.3;
obj_y(3) = 1.2;
obj_x(2) = 1;
obj_y(2) = 1;
obj_vy(1) = 2;
obj_vx(1) = 2;

% Total energy
energy = sum(obj_mass .* (obj_vx.^2 + obj_vy.^2)) / 2;

% Time step
dt = 0.01;
N_iterations = 1000;

% Create the figure
figure;

for k = 1:N_iterations
    for i = 1:N_obj
        for j = 1:N_obj
            if i ~= j
                % Calculate distance between objects
                dist = sqrt((obj_x(i) - obj_x(j))^2 + (obj_y(i) - obj_y(j))^2);
                tresh = obj_size(i) + obj_size(j);
                
                % Check if a collision occurs
                if dist <= tresh
                    % Normal vector of the collision
                    nx = (obj_x(j) - obj_x(i)) / dist;
                    ny = (obj_y(j) - obj_y(i)) / dist;
                    
                    % Relative velocity in the normal direction
                    v_rel = (obj_vx(i) - obj_vx(j)) * nx + (obj_vy(i) - obj_vy(j)) * ny;
                    
                    % Update velocities based on the conservation of momentum for elastic collisions
                    if v_rel >0  % Only update if they are moving towards each other
                        m1 = obj_mass(i);
                        m2 = obj_mass(j);
                        v1x_new = obj_vx(i) - (2 * m2 * v_rel) / (m1 + m2) * nx;
                        v1y_new = obj_vy(i) - (2 * m2 * v_rel) / (m1 + m2) * ny;
                        v2x_new = obj_vx(j) + (2 * m1 * v_rel) / (m1 + m2) * nx;
                        v2y_new = obj_vy(j) + (2 * m1 * v_rel) / (m1 + m2) * ny;
                        
                        % Assign the new velocities
                        obj_vx(i) = v1x_new;
                        obj_vy(i) = v1y_new;
                        obj_vx(j) = v2x_new;
                        obj_vy(j) = v2y_new;
                    end
                end
            end
        end
    end
    
    % Update positions based on velocities
    obj_x = obj_x + obj_vx * dt;
    obj_y = obj_y + obj_vy * dt;
    
    %wall collision
    for a=1:N_obj
    if obj_x(a)>containerSize||obj_x(a)<0
    obj_x(a)= obj_x(a) - obj_vx(a) * dt;
    obj_vx(a)=-obj_vx(a);
    end
    if obj_y(a)>containerSize||obj_y(a)<0
     obj_y(a)= obj_y(a) - obj_vy(a) * dt;
    obj_vy=-obj_vy;
    end
    end
    % Apply a condition to stop an object if its velocity is very small
    for i = 1:N_obj
        if abs(obj_vx(i)) < 1e-4 && abs(obj_vy(i)) < 1e-4
            obj_vx(i) = 0;
            obj_vy(i) = 0;
        end
    end
    
    % Clear previous frame
    clf;
    
    % Plot the objects as disks (solid circles)
    hold on;
    for i = 1:N_obj
        rectangle('Position', [obj_x(i) - obj_size(i), obj_y(i) - obj_size(i), 2*obj_size(i), 2*obj_size(i)], ...
                  'Curvature', [1, 1], 'EdgeColor', 'k', 'FaceColor', 'r');
    end
    axis([0 containerSize 0 containerSize]);
    xlabel('X Position');
    ylabel('Y Position');
    title('Collision simulation');
    daspect([1 1 1])
    grid on;
    % Pause to create animation effect
    pause(0.001); % Adjust time interval to control speed of the simulation
end
