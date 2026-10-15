%% Projecte Vehicles Task 1


%% Input data
load("model_data.mat");

%% PROCESS DATA

% 1. fixnodes for dirichlet conditions
set1=sets.subcase_0_SET_1;
set2=sets.subcase_0_SET_2;
set3=sets.subcase_0_SET_3;

sets=case_control_sets;
nnodes=size(set1,1)+size(set2,1)+size(set3,1);
fixnodes=zeros(nnodes,3);

dofDir1= getDOF(set1,3*ones(size(set1,1),1));

dir1= [set1,dofDir1,zeros(size(set1,1),1)];

dofDir21=getDOF(set2,1*ones(size(set2,1),1));
dofDir22=getDOF(set2,2*ones(size(set2,1),1));
dofDir23=getDOF(set2,6*ones(size(set2,1),1));

dir2=[set2,dofDir21,0;
      set2,dofDir22,0;
      set2,dofDir23,0];

fixnodes=[dir1;dir2];

% Index
dofs=(1:1:size(node_coords,1)*6).';
i_nd=fixnodes(:,2);
i_nn = setdiff(dofs,fixnodes(:,2));

% ud
ud=fixnodes(:,3);


% Generate Fn
g=9.81;
a_unit=[0;0;-g;0;0;0];
a=repmat(a_unit,size(node_coords,1),1);

F_aux= MAAX*a;
Fn=F_aux(i_nn);
Fdext=F_aux(i_nd);


%% 1 Computee gravity sag



function DOF = getDOF(node, localDOF)
    DOF = 6*(node - 1) + localDOF;
end