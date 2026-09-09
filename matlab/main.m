%% Clear workspace
clear; clc;
close all;

%% Variable Definition
[N_x, N_y, max_itt, rho_ref, Re, tau, maxError] = readInput("../inputs/case0.txt");

%N_x = 101;
%N_y = N_x;
%max_itt = 100;

c_s = 1/sqrt(3);
ksi = [0 1 0 -1  0 1 -1 -1  1; ...
       0 0 1  0 -1 1  1 -1 -1]';
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]';

%tau = 0.6; % Relaxation times
v = (tau - 0.5)*c_s^2; % Kinematic Viscosity
%Re = 100;
U_lid = Re*v/N_x;

%rho_ref = 2;
rho_temp = 0;

Rho = ones(1, N_y, N_x)*rho_ref;
U = zeros(2, N_y, N_x);

f = eqm_d2q9_vec(Rho, ksi, U, w);
f_new = f;
f_eq = f;

%% Main loop
%max_itt = 50000;
for itt = 1:max_itt
    for j=1:N_y
        for i=1:N_x
            if j == 1 % This is the top boundary nodes
                if i ==1 % Top-Left corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(3,j,i) = f(3,j+1,i);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(7,j,i) = f(7,j+1,i+1);
    
                    % Unknown
                    Rho_c=(Rho(1,j,i+1)+Rho(1,j+1,i))/2;
                    f_new(2,j,i) = f_new(4,j,i);
                    f_new(5,j,i) = f_new(3,j,i);
                    f_new(9,j,i) = f_new(7,j,i);
                    f_new(6,j,i) = (Rho_c - f_new(1,j,i) - f_new(3,j,i) - f_new(4,j,i) - f_new(7,j,i) - ...
                                           f_new(2,j,i) - f_new(5,j,i) - f_new(9,j,i))/2;
                    f_new(8,j,i) = f_new(6,j,i);
                elseif i == N_x % Top-right corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(3,j,i) = f(3,j+1,i);
                    f_new(6,j,i) = f(6,j+1,i-1);
    
                    % Unknown
                    Rho_c=(Rho(1,j,i-1)+Rho(1,j+1,i))/2;
                    f_new(4,j,i) = f_new(2,j,i);
                    f_new(5,j,i) = f_new(3,j,i);
                    f_new(8,j,i) = f_new(6,j,i);
                    f_new(7,j,i) = (Rho_c - f_new(1,j,i) - f_new(2,j,i) - f_new(3,j,i) - f_new(6,j,i) - ...
                                           f_new(4,j,i) - f_new(5,j,i) - f_new(8,j,i))/2;
                    f_new(9,j,i) = f_new(7,j,i);
                else % All other nodes on the top boundary
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(3,j,i) = f(3,j+1,i);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(6,j,i) = f(6,j+1,i-1);
                    f_new(7,j,i) = f(7,j+1,i+1);
    
                    % Unknown
                    f_new(5,j,i) = f_new(3,j,i);
                    Rho_b = f_new(1,j,i) + f_new(2,j,i) + f_new(4,j,i) + 2*(f_new(3,j,i) + f_new(6,j,i) + f_new(7,j,i));
                    f_new(8,j,i) = f_new(6,j,i) + (f_new(2,j,i)-f_new(4,j,i))/2 - Rho_b*U_lid/2; 
                    f_new(9,j,i) = f_new(7,j,i) - (f_new(2,j,i)-f_new(4,j,i))/2 + Rho_b*U_lid/2;
                end
            elseif j == N_y % This is the bottom boundary nodes
                if i ==1 % Bottom-Left corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(5,j,i) = f(5,j-1,i);
                    f_new(8,j,i) = f(8,j-1,i+1);
                    
                    % Unknown
                    Rho_c=(Rho(1,j,i+1)+Rho(1,j-1,i))/2;
                    f_new(2,j,i) = f_new(4,j,i);
                    f_new(3,j,i) = f_new(5,j,i);
                    f_new(6,j,i) = f_new(8,j,i);
                    f_new(7,j,i) = (Rho_c - f_new(1,j,i) - f_new(4,j,i) - f_new(5,j,i) - f_new(8,j,i) - ...
                                           f_new(2,j,i) - f_new(3,j,i) - f_new(6,j,i))/2;
                    f_new(9,j,i) = f_new(7,j,i);
                elseif i == N_x % Bottom-right corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(5,j,i) = f(5,j-1,i);
                    f_new(9,j,i) = f(9,j-1,i-1);
    
                    % Unknown
                    Rho_c=(Rho(1,j,i-1)+Rho(1,j-1,i))/2;
                    f_new(3,j,i) = f_new(5,j,i);
                    f_new(4,j,i) = f_new(2,j,i);
                    f_new(7,j,i) = f_new(9,j,i);
                    f_new(6,j,i) = (Rho_c - f_new(1,j,i) - f_new(2,j,i) - f_new(5,j,i) - f_new(9,j,i) - ...
                                           f_new(3,j,i) - f_new(4,j,i) - f_new(7,j,i))/2;
                    f_new(8,j,i) = f_new(6,j,i);
                else % All other nodes on the bottom boundary
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(5,j,i) = f(5,j-1,i);
                    f_new(8,j,i) = f(8,j-1,i+1);
                    f_new(9,j,i) = f(9,j-1,i-1);
    
                    % Unknown
                    f_new(3,j,i) = f_new(5,j,i);
                    f_new(6,j,i) = f_new(8,j,i) + (f_new(4,j,i)-f_new(2,j,i))/2;
                    f_new(7,j,i) = f_new(9,j,i) - (f_new(4,j,i)-f_new(2,j,i))/2;
                end
            elseif i == 1 % This is the left boundary nodes
                f_new(1,j,i) = f(1,j,i);
                f_new(3,j,i) = f(3,j+1,i);
                f_new(4,j,i) = f(4,j,i+1);
                f_new(5,j,i) = f(5,j-1,i);
                f_new(7,j,i) = f(7,j+1,i+1);
                f_new(8,j,i) = f(8,j-1,i+1);
    
                % Unknown
                f_new(2,j,i) = f_new(4,j,i);
                f_new(6,j,i) = f_new(8,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2;
                f_new(9,j,i) = f_new(7,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2;
            elseif i == N_x % This is the right boundary nodes
                f_new(1,j,i) = f(1,j,i);
                f_new(2,j,i) = f(2,j,i-1);
                f_new(3,j,i) = f(3,j+1,i);
                f_new(5,j,i) = f(5,j-1,i);
                f_new(6,j,i) = f(6,j+1,i-1);
                f_new(9,j,i) = f(9,j-1,i-1);
    
                % Unknown
                f_new(4,j,i) = f_new(2,j,i);
                f_new(7,j,i) = f_new(9,j,i) + (f_new(5,j,i)-f_new(3,j,i))/2;
                f_new(8,j,i) = f_new(6,j,i) - (f_new(5,j,i)-f_new(3,j,i))/2;
            else  % All interior nodes
                f_new(1,j,i) = f(1,j,i);
                f_new(2,j,i) = f(2,j,i-1);
                f_new(3,j,i) = f(3,j+1,i);
                f_new(4,j,i) = f(4,j,i+1);
                f_new(5,j,i) = f(5,j-1,i);
                f_new(6,j,i) = f(6,j+1,i-1);
                f_new(7,j,i) = f(7,j+1,i+1);
                f_new(8,j,i) = f(8,j-1,i+1);
                f_new(9,j,i) = f(9,j-1,i-1);
            end
        end
    end

    % Collision
    Rho = sum(f_new, 1);
    U = pagemtimes(ksi', f_new)./Rho;

    % Feq Calculation
    f_eq = eqm_d2q9_vec(Rho, ksi, U, w);

    % BGK + Update
    f = f_new - (f_new - f_eq)/tau;
    
    if (mod(itt, 10) == 0)
        fprintf("Itt: %i \n", itt);
    end
end

writematrix(flipud(Rho), "output/rho.dat", "Delimiter", ' ');


%data = load("Ghia_Re100.mat");

% figure
% contourf(flipud(squeeze(Rho)),3)
% axis equal tight 
% 
% figure 
% quiver(flipud(squeeze(U(1,:,:))), flipud(squeeze(U(2,:,:))), 5)
% axis equal tight 

% figure;
% hold on;
% ux = squeeze(flipud(U(1,:,51)))';
% plot(data.u_Ghia, data.y_Ghia,  "R");
% plot( ux/max(abs(ux)), flip(linspace(0, 1, 101)), "--B")
% hold off;
% 
% figure;
% hold on;
% uy = squeeze(flipud(U(2,51,:)))';
% plot(data.x_Ghia, data.v_Ghia/max(abs(data.v_Ghia)),  "R");
% plot(linspace(0, 1, 101), uy/max(abs(uy)), "--B")
% hold off;