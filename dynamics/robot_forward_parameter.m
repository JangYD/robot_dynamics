clc;
clear;
format long

%%standard dh parameter
% a = [0 -0.6127 -0.57155 0 0 0];
% al = [pi/2 0 0 pi/2 -pi/2 0];
% d = [0.1807 0 0 0.17415 0.11985 0.11655];
% th = [0 0 0 0 0 0];

%%modified dh parameter
% a = [0 0 -0.6127 -0.57155 0 0];
% al = [0 pi/2 0 0 pi/2 -pi/2 ];
% d = [0.1807 0 0 0.17415 0.11985 0.11655];
% th = [0 0 0 0 0 0];




% a = [ 0 -0.6127 -0.57155 ];
% al = [0 0 pi];
% d = [0.1807 0 0.17415];
% % th = sym('q', [1 3], 'real');
% th = [0 0 0];

a = [ 0 ];
al = [pi/2];
d = [0.1807];
th = [0 ];







dh_param = dh_parameter(a, al, d, th);

robot = make_robot(dh_param);
% showdetails(robot);
% show(robot);

% test = [     0.287459277046154   0.325235731089429  -2.323583782157943  -0.258929046128030   1.385557486771657 1.349704938089496];
% test = [     0.263003807022697   0.318225126668248  -2.083767463389982  -1.376049488970811   1.307792012600387  -1.570795984868527];
% test = [       0.263004082831239   0.318225157529930  -2.083766510664403  -1.376051300455320   1.307792243963658  -1.570796326794897 ];

 % test = [0 0 -3*pi/4 0 0 0];


% robot.Bodies{1}.Joint.HomePosition = 0;
% robot.Bodies{2}.Joint.HomePosition =  -1.815730667858608; 
% robot.Bodies{3}.Joint.HomePosition =  0;
% robot.Bodies{4}.Joint.HomePosition =  2.001415346416216   ;+.k,/
% robot.Bodies{5}.Joint.HomePosition =1.385111648237289 ;
% robot.Bodies{6}.Joint.HomePosition = 0.785398163397448;
for i = 1 : 1
    robot.Bodies{i}.Joint.HomePosition = th(i);
end

th = [robot.Bodies{1}.Joint.HomePosition %robot.Bodies{2}.Joint.HomePosition...
       % robot.Bodies{3}.Joint.HomePosition %robot.Bodies{4}.Joint.HomePosition...
      % robot.Bodies{5}.Joint.HomePosition %robot.Bodies{6}.Joint.HomePosition...
      ];
robot_dh_param = dh_parameter(a, al, d, th)

[robot_end_position, link_position] = forward_kinematics(robot_dh_param);

robot_end_matrix = robot_end_position

show(robot)

function [robot_end_position, link_position, R_matrix] = forward_kinematics_mod(dh_parameter)
    n = size(dh_parameter, 1);
    robot_end_position = eye(4);
    link_position = cell(n, 1);
    R_matrix = cell(n, 1);

    for i = 1:n
        a = dh_parameter(i,1); alpha = dh_parameter(i,2);
        d = dh_parameter(i,3); theta = dh_parameter(i,4);

        cos_theta = cos(theta); 
        sin_theta = sin(theta);
        cos_alpha = cos(alpha); 
        sin_alpha = sin(alpha);

        A = [cos_theta,             -sin_theta,            0,                   a;
             sin_theta*cos_alpha,   cos_theta*cos_alpha, -sin_alpha, -d*sin_alpha;
             sin_theta*sin_alpha,   cos_theta*sin_alpha,  cos_alpha,  d*cos_alpha;
             0,                      0,                    0,                1];
        
        robot_end_position = robot_end_position * A;
        link_position{i} = robot_end_position;
        R_matrix{i} = A;
    end
end