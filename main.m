%% Projecte Vehicles Task 1

clc; clear;

%% Load data
load("model_data.mat");
K=KAAX; % units are in N,mm,tn
M=MAAX;

%% DEFINE GENERAL BC & INPUTS

g      = 9810; % mm/s^2
nnodes = size(node_coords,1);
dofs   = (1:1:nnodes*6).';
ndofs  = size(dofs,1);

% 1. fixnodes for dirichlet conditions
sets.s1=case_control_sets.subcase_0_SET_1;
sets.s2=case_control_sets.subcase_0_SET_2;
sets.s3=case_control_sets.subcase_0_SET_3;

% Create fixed nodes matrix
dir1 = [sets.s1,getDOF(sets.s1,3*ones(size(sets.s1,1),1)), zeros(size(sets.s1,1),1)];
dir2 = [sets.s2,getDOF(sets.s2,1*ones(size(sets.s2,1),1)),0;
        sets.s2,getDOF(sets.s2,2*ones(size(sets.s2,1),1)),0;
        sets.s2,getDOF(sets.s2,6*ones(size(sets.s2,1),1)),0];

fixnodes=[dir1;dir2];

% Index
i_nd = fixnodes(:,2);
i_nn = setdiff(dofs,fixnodes(:,2));

% ud
ud = fixnodes(:,3);


%% 1A COMPUTE R FORCES AXIAL SUPPORT AND CENTRAL MEMBRANE

% Extract K submatrices
Kdd = K(i_nd, i_nd);
Kdn = K(i_nd, i_nn);
Knd = K(i_nn, i_nd);
Knn = K(i_nn, i_nn);

% Z DIRECTION
a     = repmat([0;0;-g;0;0;0],nnodes,1);
F_aux = M*a;
Fn    = F_aux(i_nn);
Fdext = F_aux(i_nd);

un   = Knn\(Fn-Knd*ud);
Fd   = Kdd*ud+Kdn*un;
gravZ.Rf = Fd-Fdext;
gravZ.R_axial_      = gravZ.Rf(1); % reaction of axial support alone
gravZ.R_axialTotal  = sum(gravZ.Rf(1:size(dir1,1))); % total reaction of the axial supports
gravZ.Rx_membrane   = gravZ.Rf(size(dir1,1)+1);
gravZ.Ry_membrane   = gravZ.Rf(size(dir1,1)+2);
gravZ.Mz_membrane   = gravZ.Rf(size(dir1,1)+3);

gravZ.U         = zeros(ndofs,1);
gravZ.U(i_nd) = ud;
gravZ.U(i_nn) = un;


% Y DIRECTION
a     = repmat([0;-g;0;0;0;0],nnodes,1);
F_aux = M*a;
Fn    = F_aux(i_nn);
Fdext = F_aux(i_nd);

un    = Knn\(Fn-Knd*ud);
Fd    = Kdd*ud+Kdn*un;
gravY.Rf  = Fd-Fdext;
gravY.R_axial_      = gravY.Rf(1); % reaction of axial support alone
gravY.R_axialTotal  = sum(gravY.Rf(1:size(dir1,1))); % total reaction of the axial supports
gravY.Rx_membrane   = gravY.Rf(size(dir1,1)+1);
gravY.Ry_membrane   = gravY.Rf(size(dir1,1)+2);
gravY.Mz_membrane   = gravY.Rf(size(dir1,1)+3);

gravY.U       = zeros(ndofs,1);
gravY.U(i_nd) = ud;
gravY.U(i_nn) = un;


% Print reactions results
fprintf('Gravity in z: \n');
fprintf('Total axial support z-reaction = %.6g N\n', gravZ.R_axialTotal);
fprintf('Total membrane x-reaction = %.6g N\n', gravZ.Rx_membrane);
fprintf('Total membrane y-reaction = %.6g N\n', gravZ.Ry_membrane);

fprintf('\nGravity in y: \n');
fprintf('Total axial support z-reaction = %.6g N\n', gravY.R_axialTotal);
fprintf('Total membrane x-reaction = %.6g N\n', gravY.Rx_membrane);
fprintf('Total membrane y-reaction = %.6g N\n', gravY.Ry_membrane);

%% 1B COMPUTE MASS AND REACTION FORCES ARE SAME IN BOTH DIRECTIONS

% Total mass with mass matrix
a_aux = repmat([1;0;0;0;0;0],nnodes,1);
mTotal= a_aux.'*M*a_aux;  % we project the mass matrix into a rb space

% print Mass results
fprintf('\nMass Comparison:\n');
fprintf('M from Mass matrix = %.6g tn\n', mTotal);
fprintf('M from Rf in z direction = %.6g tn\n', gravZ.R_axialTotal/g);
fprintf('M from Rf in y direction = %.6g tn\n', gravY.Ry_membrane/g);


%% 1C COMPUTE RMS SFE

% Optical surface nodes
optNodes = sets.s3;
x = node_coords(optNodes,1);
y = node_coords(optNodes,2);

% We extract the vertical displacements of the surface nodes
wZ = gravZ.U(getDOF(optNodes,3));
wY = gravY.U(getDOF(optNodes,3));

A = [ones(size(optNodes,1),1), x, y];

% Z DIRECTION
coeff_z     = A \ wZ;                % We fit the displacements in a plane   
wZ_fitted = A*coeff_z;               % Displacements from the fitted plane
eZ        = wZ - wZ_fitted;          % Residual
RMS_z     = sqrt(mean(eZ.^2))*1e6;

% Y DIRECTION 
coeff_y       = A \ wY;
wY_fitted = A*coeff_y;
eY        = wY - wY_fitted;
RMS_y     = sqrt(mean(eY.^2))*1e6;

fprintf('\nRMS results:\n');
fprintf('RMS in z direction = %.6g nm\n', RMS_z);
fprintf('RMS in y direction = %.6g nm\n', RMS_y);

%% 2. POLISHING ORIENTATION

% Range of angles to evaluate (from 0 to 90 degrees)
alpha_vec = 0:1:90; 
n_alpha = length(alpha_vec);


calc_rms = @(disp) sqrt(mean((disp - A*(A\disp)).^2)) * 1e6;

% Pre-allocate memory to store results
avg_rms = zeros(n_alpha, 1);     
rms_opt = zeros(n_alpha, 1);    
rms_zen = zeros(n_alpha, 1);     
rms_hor = zeros(n_alpha, 1);     

% Main loop: test every possible polishing angle (alpha_p)
for i = 1:n_alpha
    alpha_p = alpha_vec(i);

    % Mirror displacement if polished at alpha_p
    w_polish = wZ * sind(alpha_p) + wY * cosd(alpha_p);

    % Temporary array to store errors at each observation angle
    rms_obs = zeros(n_alpha, 1);

    % Inner loop: observe from all angles in the sky (alpha)
    for i = 1:n_alpha
        alpha_obs = alpha_vec(i);

        % Mirror displacement when observing at alpha_obs
        w_obs = wZ * sind(alpha_obs) + wY * cosd(alpha_obs);

        % SFE = Difference in shape (removing plane with calc_rms)
        rms_obs(i) = calc_rms(w_obs - w_polish);
    end

    % Save the average error for this polishing angle
    avg_rms(i) = mean(rms_obs);

    % Save specific curves requested in part 
    if alpha_p == 90
        rms_zen = rms_obs;
    elseif alpha_p == 0
        rms_hor = rms_obs;
    end
end

% Find the optimal polishing angle 
[min_avg_val, opt_idx] = min(avg_rms);
opt_alpha_p = alpha_vec(opt_idx);

% Recalculate the error curve for the optimal polishing 
w_polish_opt = wZ * sind(opt_alpha_p) + wY * cosd(opt_alpha_p);
for i = 1:n_alpha
    w_obs = wZ * sind(alpha_vec(i)) + wY * cosd(alpha_vec(i));
    rms_opt(i) = calc_rms(w_obs - w_polish_opt);
end
max_rms_opt = max(rms_opt);


fprintf('Optimal polishing elevation (alpha_p) = %d deg\n', opt_alpha_p);
fprintf('Average SFE for optimal polish = %.2f nm-rms\n', min_avg_val);
fprintf('Maximum SFE for optimal polish = %.2f nm-rms\n', max_rms_opt);
fprintf('\nComparison:\n');
fprintf('Polish at Zenith (90 deg): Avg SFE = %.2f nm\n', avg_rms(end));
fprintf('Polish at Horizon (0 deg): Avg SFE = %.2f nm\n', avg_rms(1));

% Plot the RMS SFE versus alpha
figure('Name', 'SFE vs Elevation');
plot(alpha_vec, rms_opt, 'g-', 'LineWidth', 2); hold on;
plot(alpha_vec, rms_zen, 'b--', 'LineWidth', 1.5);
plot(alpha_vec, rms_hor, 'r:', 'LineWidth', 1.5);
xlabel('Telescope elevation \alpha (degrees)');
ylabel('RMS SFE (nm)');
title('Surface Figure Error vs Elevation');
legend(sprintf('Optimal Polish (\\alpha_p = %d^\\circ)', opt_alpha_p), ...
    'Polish at Zenith (90^\\circ)', ...
    'Polish at Horizon (0^\\circ)', 'Location', 'best');
grid on;


%% 3. Eigen modes

% RB modes unconstrained
[V, D] = eigs(K, M, 7, 'smallestabs'); 
freq = sqrt(D)./(2*pi);
freq = real(diag(freq));
[freq,~] = sort(freq);

% RB modes constrained with BC ex 1
Mnn=M(i_nn, i_nn);
[Phi, D_c] = eigs(Knn, Mnn, 7, 'smallestabs'); 

mode6 = zeros(ndofs,1);
mode6(i_nn) = Phi(:, 6);
Vplot=reshape(mode6,[],6);
fillhdf('template.h5','6_eigenmode.h5',Vplot);

lambda = real(diag(D_c));
[lambda,idx] = sort(lambda);

Freq_dyn = sqrt(lambda)/(2*pi);
wn = sqrt(lambda);


%% 4. Dynamic SFE tip-tilt motion

% Construction of acceleration vector for sinusoidal movement

a = zeros(ndofs, 6);
a(1:6:end, 1) = 1;
a(2:6:end, 2) = 1;
a(3:6:end, 3) = 1;
% We solve the eigenmode problem
[Phi, D_c] = eigs(Knn, Mnn, 7, 'smallestabs'); 

theta0=100e-6; %rad
zeta = 0.02;    % amortiguament 2%
f=1e-03:500;
w= 2*pi*f;

% Calculate K and M for each mode
mModal = zeros(size(Phi,2),1);
kModal = zeros(size(Phi,2),1);
b = zeros(size(Phi,2),a1);

for j = 1:size(Phi,2)
    phi = Phi(:,j);
    mModal(j) = phi'*Mnn*phi;
    kModal(j) = phi'*Knn*phi;
end
Fxi=Phi'*Fn;

% Damping modal freq response
Xi=zeros(7,length(w));
for n=1:length(w)
    alpha0 = theta0*w(n)^2;  % amplitut angular acceleration
    % 2. Rigid-body acceleration
    a_rb = zeros(ndofs,1);
    a_rb(2:6:end) = -alpha0*node_coords(:,3);
    a_rb(3:6:end) =  alpha0*node_coords(:,2);
    a_rb(4:6:end) =  alpha0;
    Xi(:,n)= Fxi./(-w(n)^2.*mModal+kModal+1i.*2*mModal.*wn.*zeta.*w(n));
end