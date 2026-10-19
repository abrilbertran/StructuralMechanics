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
z.Rf = Fd-Fdext;
z.R_axial_      = z.Rf(1); % reaction of axial support alone
z.R_axialTotal  = sum(z.Rf(1:size(dir1,1))); % total reaction of the axial supports
z.Rx_membrane   = z.Rf(size(dir1,1)+1);
z.Ry_membrane   = z.Rf(size(dir1,1)+2);
z.Mz_membrane   = z.Rf(size(dir1,1)+3);

% Y DIRECTION
a     = repmat([0;-g;0;0;0;0],nnodes,1);
F_aux = M*a;
Fn    = F_aux(i_nn);
Fdext = F_aux(i_nd);

un    = Knn\(Fn-Knd*ud);
Fd    = Kdd*ud+Kdn*un;
y.Rf  = Fd-Fdext;
y.R_axial_      = y.Rf(1); % reaction of axial support alone
y.R_axialTotal  = sum(y.Rf(1:size(dir1,1))); % total reaction of the axial supports
y.Rx_membrane   = y.Rf(size(dir1,1)+1);
y.Ry_membrane   = y.Rf(size(dir1,1)+2);
y.Mz_membrane   = y.Rf(size(dir1,1)+3);

% Print reactions results
fprintf('Gravity in z: \n');
fprintf('Total axial support z-reaction = %.6g N\n', z.R_axialTotal);
fprintf('Total membrane x-reaction = %.6g N\n', z.Rx_membrane);
fprintf('Total membrane y-reaction = %.6g N\n', z.Ry_membrane);

fprintf('\nGravity in y: \n');
fprintf('Total axial support z-reaction = %.6g N\n', y.R_axialTotal);
fprintf('Total membrane x-reaction = %.6g N\n', y.Rx_membrane);
fprintf('Total membrane y-reaction = %.6g N\n', y.Ry_membrane);

%% 1B COMPUTE MASS AND REACTION FORCES ARE SAME IN BOTH DIRECTIONS

% Total mass with mass matrix
a_aux = repmat([1;0;0;0;0;0],nnodes,1);
mTotal= a_aux.'*M*a_aux;  % we project the mass matrix into a rb space

% print Mass results
fprintf('\nMass Comparison:\n');
fprintf('M from Mass matrix = %.6g tn\n', mTotal);
fprintf('M from Rf in z direction = %.6g tn\n', z.R_axialTotal/g);
fprintf('M from Rf in y direction = %.6g tn\n', y.Ry_membrane/g);


%% 1C COMPUTE RMS SFE

