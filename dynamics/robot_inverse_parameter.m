clc;
clear;

syms q1 q2 q3 q4 q5 q6 real;
a = [0 -0.6127 -0.57155 0 0 0];
al = [pi/2 0 0 pi/2 -pi/2 0];
d = [0.1807 0 0 0.17415 0.11985 0.11655];
th = [q1 q2 q3 q4 q5 q6];

dh_param = dh_parameter(a, al, d, th);


robot_end_mat = [      
  -0.707106781186547   0.707106781186548                   0  -0.993299866938037
  -0.000000000000000  -0.000000000000000  -1.000000000000000  0.209300000000000
  -0.707106781186548  -0.707106781186547   0.000000000000000   0.969593628512379
                   0                   0                   0   1.000000000000000
];

q = inverse_kinematics(dh_param, robot_end_mat)
theta = zeros(1,5);
for i = 1 : 6
    theta(i) = q(i,8);
end


inv_dh_param = dh_parameter(a, al, d, theta);

[robot_end_position, link_position] = forward_kinematics(inv_dh_param);
robot_mat = robot_end_position
robot = make_robot(inv_dh_param);
robot.Bodies{1}.Joint.HomePosition = theta(1);
robot.Bodies{2}.Joint.HomePosition = theta(2);
robot.Bodies{3}.Joint.HomePosition = theta(3);
robot.Bodies{4}.Joint.HomePosition = theta(4);
robot.Bodies{5}.Joint.HomePosition = theta(5);
robot.Bodies{6}.Joint.HomePosition = theta(6);
show(robot);




