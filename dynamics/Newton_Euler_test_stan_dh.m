clc;
clear;
digitsOld = digits(10);

%%
q_dot = sym('q_dot_',[1,6], 'real');
q_dot_dot = sym('q_dot_dot_',[1,6], 'real');
                                
a = [0 -0.6127 -0.57155 0 0 0];
al = [pi/2 0 0 pi/2 -pi/2 0];
d = [0.1807 0 0 0.17415 0.11985 0.11655];
th = sym('q', [1 6], 'real');

mass = [7.369 13.051 3.989 2.1 1.98 0.615];
g = [0; 0; 9.81];

robot_dh_parameter = dh_parameter(a, al, d, th);
[robot_end_matrix, robot_link_matrix, R_matrix] = forward_kinematics(robot_dh_parameter);

robot_length = length(robot_link_matrix);

center_of_mass = cell(robot_length,1);
center_of_mass{1} = [0.021; 0.000; 0.027];
center_of_mass{2} = [0.38; 0.000; 0.158];
center_of_mass{3} = [0.24; 0.000; 0.068];
center_of_mass{4} = [0.000; 0.007; 0.018];
center_of_mass{5} = [0.000; 0.007; 0.018];
center_of_mass{6} = [0; 0; -0.026];

inertia_matrix = cell(robot_length,1);
inertia_matrix{1} = [0.0341, 0.0000, -0.0043; 0.0000, 0.0353, 0.0001; -0.0043, 0.0001, 0.0216];
inertia_matrix{2} = [0.0281, 0.0001, -0.0156; 0.0001, 0.7707, 0.0000; -0.0156, 0.0000, 0.7694];
inertia_matrix{3} = [0.0101, 0.0001, 0.0092; 0.0001, 0.3093, 0.0000; 0.0092, 0.0000, 0.3065];
inertia_matrix{4} = [0.0030, -0.0000, -0.0000; -0.0000, 0.0022, -0.0002; -0.0000, -0.0002, 0.0026];
inertia_matrix{5} = [0.0030, -0.0000, -0.0000; -0.0000, 0.0022, -0.0002; -0.0000, -0.0002, 0.0026];
inertia_matrix{6} = [0.0000, 0.0000, -0.0000; 0.0000, 0.0004, 0.0000; -0.0000, 0.0000, 0.0003];

H_mat = robot_link_matrix;

R{7} = [1 0 0; 0 1 0; 0 0 1];

rot_x = cell(robot_length, 1);
for i = 1 : robot_length
    H_rot_x = Rx(al(i));
    rot_x{i} = round(H_rot_x(1:3,1:3));
end

H_rot_x = cell(robot_length, 1);
for i = 1 : robot_length
    H_rot_x{i} = eye(3);
    for j = 1 : i
        H_rot_x{i} = H_rot_x{i}*rot_x{j};
    end
end

link = cell(robot_length,1);
for i = 1 : robot_length
    link{i} = rot_x{i}'*[a(i); 0; d(i);];
end

IG = cell(robot_length,1);
for i = 1 : robot_length
    R{i} = R_matrix{i}(1:3,1:3);
    IG{i} = rot_x{i}*inertia_matrix{i}*rot_x{i}';
end

Z = cell(robot_length,1);
z0 = [0;0;1];
for i = 1 :robot_length
     Z{i} =  rot_x{i}'*z0;
end

%%
com1_front = [0.021;0.1807;0.027];
com1_back = [0.021;0;0.027];

com2_front = [-0.2327;0;0.158];
com2_back = [0.38;0;0.158];
 
com3_front = [-0.33155; 0; 0.068];
com3_back = [0.24; 0; 0.068];

com4_front = [0.000; 0.18115; 0.018];
com4_back = [0.000; 0.007; 0.018];

com5_front = [0;  -0.11285; 0.018];
com5_back = [0.000; 0.007; 0.018];

com6_front =[0; 0; 0.09055];
com6_back = [0; 0; -0.026];

w0 = [0;0;0];
w1 = R{1}'*(w0 + q_dot(1)*z0);
w2 = R{2}'*(w1 + q_dot(2)*z0);
w3 = R{3}'*(w2 + q_dot(3)*z0);
w4 = R{4}'*(w3 + q_dot(4)*z0);
w5 = R{5}'*(w4 + q_dot(5)*z0);
w6 = R{6}'*(w5 + q_dot(6)*z0);

acc_w0 = [0;0;0];
acc_w1 = R{1}'*acc_w0 + cross(w1, R{1}'*q_dot(1)*z0) + R{1}'*q_dot_dot(1)*z0;
acc_w2 = R{2}'*acc_w1 + cross(w2, R{2}'*q_dot(2)*z0) + R{2}'*q_dot_dot(2)*z0;
acc_w3 = R{3}'*acc_w2 + cross(w3, R{3}'*q_dot(3)*z0) + R{3}'*q_dot_dot(3)*z0;
acc_w4 = R{4}'*acc_w3 + cross(w4, R{4}'*q_dot(4)*z0) + R{4}'*q_dot_dot(4)*z0;
acc_w5 = R{5}'*acc_w4 + cross(w5, R{5}'*q_dot(5)*z0) + R{5}'*q_dot_dot(5)*z0;
acc_w6 = R{6}'*acc_w5 + cross(w6, R{6}'*q_dot(6)*z0) + R{6}'*q_dot_dot(6)*z0;

acc_v0 = [0;0;0];
acc_v1 = R{1}' * ( acc_v0 + cross(acc_w0, [0;0;0]) + cross(w0, cross(w0, [0;0;0])) );
acc_v2 = R{2}' * ( acc_v1 + cross(acc_w1, link{1}) + cross(w1, cross(w1, link{1})) );
acc_v3 = R{3}' * ( acc_v2 + cross(acc_w2, link{2}) + cross(w2, cross(w2, link{2})) );
acc_v4 = R{4}' * ( acc_v3 + cross(acc_w3, link{3}) + cross(w3, cross(w3, link{3})) );
acc_v5 = R{5}' * ( acc_v4 + cross(acc_w4, link{4}) + cross(w4, cross(w4, link{4})) );
acc_v6 = R{6}' * ( acc_v5 + cross(acc_w5, link{5}) + cross(w5, cross(w5, link{5})) );


acc_vc1 = acc_v1 + cross(acc_w1,com1_front) + cross(w1, cross(w1,com1_front));
acc_vc2 = acc_v2 + cross(acc_w2,com2_front) + cross(w2, cross(w2,com2_front));
acc_vc3 = acc_v3 + cross(acc_w3,com3_front) + cross(w3, cross(w3,com3_front));
acc_vc4 = acc_v4 + cross(acc_w4,com4_front) + cross(w4, cross(w4,com4_front));
acc_vc5 = acc_v5 + cross(acc_w5,com5_front) + cross(w5, cross(w5,com5_front));
acc_vc6 = acc_v6 + cross(acc_w6,com6_front) + cross(w6, cross(w6,com6_front));

g1 = H_mat{1}(1:3,1:3)'*g;
g2 = H_mat{2}(1:3,1:3)'*g;
g3 = H_mat{3}(1:3,1:3)'*g;
g4 = H_mat{4}(1:3,1:3)'*g;
g5 = H_mat{5}(1:3,1:3)'*g;
g6 = H_mat{6}(1:3,1:3)'*g;

f7 = [0;0;0];
f6 = R{7}*f7 + mass(6)*acc_vc6 - mass(6)*g6;
f5 = R{6}*f6 + mass(5)*acc_vc5 - mass(5)*g5;
f4 = R{5}*f5 + mass(4)*acc_vc4 - mass(4)*g4;
f3 = R{4}*f4 + mass(3)*acc_vc3 - mass(3)*g3;
f2 = R{3}*f3 + mass(2)*acc_vc2 - mass(2)*g2;
f1 = R{2}*f2 + mass(1)*acc_vc1 - mass(1)*g1;

n7 = [0;0;0];
n6 = IG{6}*acc_w6 + cross(w6, IG{6}*w6) - cross(f6, com6_front) + cross(R{7}*f7, com6_back) + R{7}*n7;
n5 = IG{5}*acc_w5 + cross(w5, IG{5}*w5) - cross(f5, com5_front) + cross(R{6}*f6, com5_back) + R{6}*n6;
n4 = IG{4}*acc_w4 + cross(w4, IG{4}*w4) - cross(f4, com4_front) + cross(R{5}*f5, com4_back) + R{5}*n5;
n3 = IG{3}*acc_w3 + cross(w3, IG{3}*w3) - cross(f3, com3_front) + cross(R{4}*f4, com3_back) + R{4}*n4;
n2 = IG{2}*acc_w2 + cross(w2, IG{2}*w2) - cross(f2, com2_front) + cross(R{3}*f3, com2_back) + R{3}*n3;
n1 = IG{1}*acc_w1 + cross(w1, IG{1}*w1) - cross(f1, com1_front) + cross(R{2}*f2, com1_back) + R{2}*n2;

tau1 = n1'*(R{1}'*z0);
tau2 = n2'*(R{2}'*z0);
tau3 = n3'*(R{3}'*z0);
tau4 = n4'*(R{4}'*z0);
tau5 = n5'*(R{5}'*z0);
tau6 = n6'*(R{6}'*z0);

t1 = vpa(simplify(subs(tau1, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
t2 = vpa(simplify(subs(tau2, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
t3 = vpa(simplify(subs(tau3, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
t4 = vpa(simplify(subs(tau4, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
t5 = vpa(simplify(subs(tau5, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))
t6 = vpa(simplify(subs(tau6, [th(1) th(2) th(3) th(4) th(5) th(6)], [0 pi/4 0 2*pi/5 5*pi/6 0])))

t1 = simplify(t1)
t2 = simplify(t2)
t3 = simplify(t3)
t4 = simplify(t4)
t5 = simplify(t5)
t6 = simplify(t6)
