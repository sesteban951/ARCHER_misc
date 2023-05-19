%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         MODIFY CLF and CBF values in YAML config
%%%   note: this script does not preserve comments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

% add YAML repo path
repo_path = "../yaml";
addpath(genpath('/yaml'))

% Linear System after Feedback linearization
eta_dim = 6;
d = eta_dim/2;

F = [zeros(d,d), eye(d,d);
     zeros(d,d), zeros(d,d)];

G = [zeros(d,d);
     eye(d,d)];

% controllabilty (always controllable)
cont = ctrb(F,G);
controllable = (eta_dim == rank(cont))

% select controller design
control_design = 1;   % 1 - LQR gains
                      % 2 - pole placements
                      % 3 - custom gains

if control_design == 1
    % LQR gains
    q_xi = 1;
    q_w = 1;
    Qx = diag([q_xi q_xi q_xi q_w q_w q_w]);
    
    qu = 1;
    Qu = diag([qu qu qu]);
    
    K_lqr = lqr(F,G,Qx,Qu)
    Acl = F-G*K_lqr
    Acl_eig = eig(Acl)
    A_hurwitz = all(Acl_eig < 0)

elseif control_design == 2
    % pole placemnt
    pol = [-1 + 0i;
           -1 + 0i;
           -1 + 0i;
           -2 + 0i;
           -2 + 0i;
           -2 + 0i];
    K_pl = place(F,G,pol)
    Acl = F-G*K_pl
    Acl_eig = eig(Acl)
    A_hurwtiz = all(Acl_eig < 0)

elseif control_design == 3
    % arbitrary Gains
    kp = 1;
    Kp = diag([kp kp kp]);
    kd = 1;
    Kd = diag([kd kd kd]);
    K = [Kp Kd]
    
    Acl = F-G*K
    Acl_eig = eig(Acl)
    A_hurwitz = all(Acl_eig < 0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%% THINGS TO ADD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lyapunov rate
lam = 1;

% CTLE, dealing with 6x6 matrices
Q = diag([1 1 1 1 1 1])
P = lyap(Acl', Q)

% CBF params
alph1 = 1.0;
alph2 = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% YAML stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLF struct
field1 = 'Q';  value1 = vec(Q);
field2 = 'P';  value2 = vec(P);
field3 = 'lamd';  value3 = lam;

s_CLF = struct(field1, value1,...
               field2, value2,...
               field3, value3);

% CBF struct
field4 = 'alph1'; value4 = alph1;
field5 = 'alph2'; value5 = alph2;
field6 = 'null'; value6 = [0 0 0]; % so you dont get weird yaml formatting

s_CBF = struct(field4, value4,...
               field5, value5,...
               field6, value6);

% get old YAML file data
data = yaml.loadFile("gains.yaml"); % loads in struct full of structs
data.CLF_CBF.CLF = s_CLF;  % add CLF data
data.CLF_CBF.CBF = s_CBF;  % add CBF data

% write to both old and new data into new file (deletes comments)
yaml.dumpFile("gains.yaml",data)

% check by loading and reading file
result = yaml.loadFile("gains.yaml")
result = yaml.dump(result)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vectorize matrix (columnwise stack)
function v = vec(M)

    v = [];
    [~, col] = size(M);

    for j = 1:col
        v = [v; M(:,j)];
    end

end






