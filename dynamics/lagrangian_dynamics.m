 clc;
clear;
format long

% digitsOld = digits(14);

q_dot = sym('q_dot_',[1,6], 'real');
q_dot_dot = sym('q_dot_dot_',[1,6], 'real');

a = [0 -0.6127 -0.57155 0 0 0];
al = [pi/2 0 0 pi/2 -pi/2 0];
d = [0.1807 0 0 0.17415 0.11985 0.11655];
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
robot_center_of_mass{1} = [0.021; 0.1807; 0.027];
robot_center_of_mass{2} = [-0.2327; 0.000; 0.158];
robot_center_of_mass{3} = [-0.33155; 0.000; 0.068];
robot_center_of_mass{4} = [0.000; 0.18115; 0.018];
robot_center_of_mass{5} = [0.000; -0.11285; 0.018];
robot_center_of_mass{6} = [0; 0; 0.09055];

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

H_mat = robot_link_matrix;

R = cell(robot_length, 1);
T = cell(robot_length, 1);
for i = 1 :robot_length
    R{i} = robot_link_matrix{i}(1:3,1:3);
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
    Jacob_com_pos{i} = (H_mat{i}(1:3,1:3)*robot_center_of_mass{i}) +link_pos ;

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
    Ke{i} = Ke{i} + mass(i)*(Jvcom{i}'*Jvcom{i}) + Jw{i}'*(H_mat{i}(1:3,1:3)*inertia_matrix{i}*H_mat{i}(1:3,1:3)')*Jw{i};
end
vpa(simplify(Jvcom{1}'*Jvcom{1}))
vpa(simplify(Jw{1}'*(H_mat{1}(1:3,1:3)*inertia_matrix{1}*H_mat{1}(1:3,1:3)')*Jw{1}))
vpa(simplify(Jvcom{2}'*Jvcom{2}))
vpa(simplify(Jw{2}'*(H_mat{2}(1:3,1:3)*inertia_matrix{2}*H_mat{2}(1:3,1:3)')*Jw{2}))

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
    for k = 1 : n
        for j = 1 : n
            
            for i = 1 : n
                C{n}(k,j) = C{n}(k,j)+ 0.5*(diff(M{n}(k,j), th(i)) + diff(M{n}(k,i), th(j)) - diff(M{n}(i,j), th(k)))*q_dot(i);
            end
        end
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

% M6_subs = vpa(simplify(subs(M{6}, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
% C6_subs = vpa(simplify(subs(C{6}, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
% G6_subs = vpa(simplify(subs(G{6}, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))

