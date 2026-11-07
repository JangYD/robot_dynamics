clc;
clear;
format long

q_dot = sym('q_dot_',[1,6], 'real');
q_dot_dot = sym('q_dot_dot_',[1,6], 'real');

a = [  0 -0.6127 -0.57155 0 0 0];
al = [  pi/2 0 0  pi/2  -pi/2 0];
d = [ 0.1807 0 0  0.17415  0.11985  0.11655];
th = sym('q', [1 6], 'real');


robot_dh_parameter = dh_parameter(a, al, d, th);
[robot_end_matrix, robot_link_matrix, R_matrix] = forward_kinematics(robot_dh_parameter);
robot_length = length(robot_link_matrix);

mass = sym('mass', [1 6], 'real');
mass(1) = 7.369;
mass(2) = 13.051; 
mass(3) = 3.989;
mass(4) = 2.1;
mass(5) = 1.98;
mass(6) = 0.615;

g = [0; 0; -9.81];

robot_center_of_mass = cell(robot_length, 1);
robot_center_of_mass{1} = [0.021; 0; 0.027];
robot_center_of_mass{2} = [0.38; 0; 0.158];
robot_center_of_mass{3} = [0.24; 0; 0.068];
robot_center_of_mass{4} = [0; 0.007; 0.018];
robot_center_of_mass{5} = [0; 0.007; 0.018];
robot_center_of_mass{6} = [0; 0; -0.026];

inertia_matrix = cell(robot_length,1);
inertia_matrix{1} = [0.0341, 0.0000, -0.0043; 0.0000, 0.0353, 0.0001; -0.0043, 0.0001, 0.0216];
inertia_matrix{2} = [0.0281, 0.0001, -0.0156; 0.0001, 0.7707, 0.0000; -0.0156, 0.0000, 0.7694];
inertia_matrix{3} = [0.0101, 0.0001, 0.0092; 0.0001, 0.3093, 0.0000; 0.0092, 0.0000, 0.3065];
inertia_matrix{4} = [0.0030, -0.0000, -0.0000; -0.0000, 0.0022, -0.0002; -0.0000, -0.0002, 0.0026];
inertia_matrix{5} = [0.0030, -0.0000, -0.0000; -0.0000, 0.0022, -0.0002; -0.0000, -0.0002, 0.0026];
inertia_matrix{6} = [0.0000, 0.0000, -0.0000; 0.0000, 0.0004, 0.0000; -0.0000, 0.0000, 0.0003];

link = cell(robot_length, 1);
for i = 1 : robot_length
    Rotx = Rx(al(i));
    link{i} = [a(i); 0; d(i);];
end

rot_x = cell(robot_length+1, 1);
for i = 1 : robot_length
    H_rot_x = Rx(al(i));
    rot_x{i} = round(H_rot_x(1:3,1:3));
end

H_mat = robot_link_matrix;

H_R_matrix = cell(robot_length, 1);
T = cell(robot_length, 1);
for i = 1 :robot_length
    H_R_matrix{i} = robot_link_matrix{i}(1:3,1:3)*rot_x{i}';
    T{i} = robot_link_matrix{i}(1:3,4);
end

Jw = cell(robot_length, 1);
Jw{1} = [0 0 1]' ;
for i = 2 : robot_length
    z = [0 0 1]';
    for j = 1 : i-1
        z = [z H_mat{j}(1:3,3)];
    end

    Jw{i} = z;
end

%% jacobian center of gravity
total_mass = 0;
Jvcom = cell(robot_length, 1);
Jacob_com_pos = cell(robot_length, 1);
robot_link_pos = cell(robot_length,1);
Jvcom_diff = cell(robot_length, 1);
Jvcom_mass = cell(robot_length, 1);

for i = 1 : robot_length
    robot_link_pos{i} = H_mat{i}(1:3,4);
end
link_pos = [0;0;0];
for i = 1 : robot_length
    total_mass = total_mass + mass(i);
    Jacob_com_pos{i} = (H_mat{i}(1:3,1:3)*robot_center_of_mass{i}) + H_mat{i}(1:3,4);

    for j = 1 : i
        Jv_com_pos_x = Jacob_com_pos{i}(1,1);
        Jv_com_pos_y = Jacob_com_pos{i}(2,1);
        Jv_com_pos_z = Jacob_com_pos{i}(3,1);
        Jvcom_diff{i} = [Jvcom_diff{i} [diff(Jv_com_pos_x, th(j)); diff(Jv_com_pos_y, th(j)); diff(Jv_com_pos_z, th(j));]];
    end

    Jvcom_diff{i} = [Jvcom_diff{i} repmat([0;0;0],1,robot_length-i)];
    Jvcom{i} = [Jvcom{i} repmat([0;0;0],1,robot_length)];
    Jvcom_mass{i} = [Jvcom_mass{i} repmat([0;0;0],1,robot_length)];

    for k = 1 : i  
        Jvcom{i} = Jvcom_diff{i};
        Jvcom_mass{i} = Jvcom_mass{i} + Jvcom_diff{k}*(mass(k)/total_mass);
    end
    link_pos = robot_link_pos{i};
end 

for i = 1 : robot_length
    Jvcom{i}(:,i+1:robot_length) = [];
    Jvcom_mass{i}(:,i+1:robot_length) = [];
end

%% Inertial
Ke = cell(robot_length,1);

for i = 1 : robot_length
    Ke{i} = [Ke{i} zeros(i,i)];
    Ke{i} = Ke{i} + mass(i)*(Jvcom{i}'*Jvcom{i}) + Jw{i}'*(H_R_matrix{i}(1:3,1:3)*inertia_matrix{i}*H_R_matrix{i}(1:3,1:3)')*Jw{i};
end

for i = 2 : robot_length
    Ke{i-1}(robot_length,robot_length) = 0;
end

M = cell(robot_length,1);
tmp = 0;
for i = 1 : robot_length
    for j = 1 : i
        tmp = tmp + Ke{j};
    end
    M{i} = tmp;
    tmp = 0;
end

%% Coriolis
C= cell(1,robot_length);
for n = 1 : robot_length
    C{n} = sym(zeros(n, n));
    Mn = M{n}(1:n,1:n);
    dM = sym(zeros(n,n,n));
    for k = 1 : n
        dM(:,:,k) = diff(Mn,th(k));
    end

    for i = 1 : n
        M_theta_i = sym(zeros(n,1));
        M_theta_j = sym(zeros(n,1));
        M_theta_k = sym(zeros(n,1));
        for j = 1 : n
            M_theta_i = M_theta_i + dM(:,i,j)*q_dot(j);
        end

        for j = 1 : n
            M_theta_j = M_theta_j + dM(:,j,i)*q_dot(j);
        end

        for j = 1 : n
            for k = 1 : n
                M_theta_k(j) = M_theta_k(j) + dM(k,i,j)*q_dot(k);
            end
        end
        C{n}(:,i) = 0.5*(M_theta_i + M_theta_j - M_theta_k);
    end
end

for i = 1 : robot_length-1
    C{i}(robot_length, robot_length) = 0;
end

%% Gravity
Pe = cell(1,robot_length);

for i = 1 : robot_length
    gravity_pos = mass(i)*g'*Jacob_com_pos{i};
    for j = 1 : i
        Pe{i} = [Pe{i} diff(gravity_pos,th(j))];
    end
end

for i = 2 : robot_length
    Pe{i-1}(1,robot_length) = 0;
end

G = cell(robot_length,1);
tmp = 0;
for i = 1 : robot_length
    for j = 1 : i
        tmp = tmp + Pe{j};
    end
    G{i} = tmp;
    tmp = 0;
end

M6_subs = vpa(simplify(subs(M{6}, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
C6_subs = vpa(simplify(subs(C{6}, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
G6_subs = vpa(simplify(subs(G{6}, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
M_last= M6_subs*q_dot_dot'
C_last = C6_subs*q_dot'
G_last = G6_subs'
tau = M_last + C_last + G_last
tau1 = simplify(tau(1))
tau2 = simplify(tau(2))
tau3 = simplify(tau(3))
tau4 = simplify(tau(4))
tau5 = simplify(tau(5))
tau6 = simplify(tau(6))

