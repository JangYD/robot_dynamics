function q = inverse_kinematics(dh_parameter, robot_end_mat) 
    %%%%%%%%%%%%%%%%%%%%%% Analytical inverse_kinematics %%%%%%%%%%%%%%%%%%%%%%%%
    a = double(dh_parameter(:,1)); 
    alpha = double(dh_parameter(:,2));
    d = double(dh_parameter(:,3));
    link = a+d;

    %%%%%%%%%%%%%%%%%%%%%% theta 1 %%%%%%%%%%%%%%%%%%%%%%%%
    P0_5 = robot_end_mat*[0; 0; -link(6); 1] - [0; 0; 0; 1];
    psi = atan2(P0_5(2,1), P0_5(1,1));
    phi = acos(link(4)/(sqrt(P0_5(1)^2 + P0_5(2)^2)));
    theta1_1 = real(psi + phi + pi/2);
    theta1_2 = real(psi - phi + pi/2);
    
    q(1,1:4) = mod(theta1_1 + pi, 2*pi) - pi;
    q(1,5:8) = mod(theta1_2 + pi, 2*pi) - pi;    

    %%%%%%%%%%%%%%%%%%%%%% theta 5 %%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : 2 
        start_ind = 4*(i-1)+1;
        H_mat0_1 =  H_mat(a(1), alpha(1), d(1), q(1, start_ind));
        H_mat1_6 = H_mat0_1\robot_end_mat;
        cos_velue = (H_mat1_6(3,4) - link(4))/link(6);
        velue = real(acos(cos_velue));
        theta5 = [velue -velue];
        q(5,4*(i-1)+1:4*(i-1)+2) = mod(theta5(1) + pi, 2*pi) - pi;
        q(5,4*(i-1)+3:4*(i-1)+4) = mod(theta5(2) + pi, 2*pi) - pi;
         
        %%%%%%%%%%%%%%%%%%%%%% theta 6 %%%%%%%%%%%%%%%%%%%%%%%%

        H_mat6_1 = inv(H_mat1_6);
        for j = 1 : 2
            start_ind = 4*(i-1) +(j-1)*j + 1;
            end_ind = start_ind +1;
            % if H_mat6_1(2,3) == 0 || H_mat6_1(1,3) == 0
            %     theta6 = real(atan2(-H_mat6_1(2,2)/sin(q(5,start_ind)), H_mat6_1(1,2)/sin(q(5,start_ind))));
            % else
            %     theta6 = real(atan2(-H_mat6_1(2,3)/sin(q(5,start_ind)), H_mat6_1(1,3)/sin(q(5,start_ind))));
            % end
            theta6 = atan2(-H_mat6_1(2,3)/sin(q(5,start_ind)), H_mat6_1(1,3)/sin(q(5,start_ind)));
            q(6,start_ind:end_ind) = mod(theta6 + pi, 2*pi) - pi;
        end
    end
    % P0_6_x = robot_end_mat(1,4);
    % P0_6_y = robot_end_mat(2,4);
    % 
    % for i = 1 : 2
    %     theta1 = q(1, 4*(i-1)+1);
    %     P1_6_z = P0_6_x*sin(theta1) - P0_6_y*cos(theta1);
    %     acos_theta5 = (P1_6_z - link(4))/link(6)
    %     % value = acos((P1_6_z - link(4))/link(6));
    % 
    %     theta5 = [value -value];
    %     q(5,4*(i-1)+1:4*(i-1)+2) = mod(theta5(1) + pi, 2*pi) - pi;
    %     q(5,4*(i-1)+3:4*(i-1)+4) = mod(theta5(2) + pi, 2*pi) - pi;
    % end
   
    %%%%%%%%%%%%%%%%%%%%%% theta 6 %%%%%%%%%%%%%%%%%%%%%%%%
    % for i = 1 : 2
    %     theta1 = q(1, 4*(i-1)+1);
    %     H_mat0_1 = H_mat(a(1), alpha(1), d(1), theta1);
    %     H_mat1_6 = H_mat0_1\robot_end_mat;
    %     H_mat6_1 = inv(H_mat1_6);
    %     for j = 1 : 2
    %         start_ind = 4*(i-1) +(j-1)*j + 1;
    %         end_ind = start_ind +1;
    %         if H_mat6_1(2,3) == 0 || H_mat6_1(1,3) == 0
    %             theta6 = real(atan2(-H_mat6_1(2,2)/sin(q(5,start_ind)), H_mat6_1(1,2)/sin(q(5,start_ind))));
    %         else
    %             theta6 = real(atan2(-H_mat6_1(2,3)/sin(q(5,start_ind)), H_mat6_1(1,3)/sin(q(5,start_ind))));
    %         end
    %         % theta6 = atan2(-H_mat6_1(2,3)/sin(q(5,start_ind)), H_mat6_1(1,3)/sin(q(5,start_ind)));
    %         q(6,start_ind:end_ind) = mod(theta6 + pi, 2*pi) - pi;
    %     end
    % end
        
    %%%%%%%%%%%%%%%%%%%%%% theta 3 %%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1 : 4
        start_ind = i*2-1;
        end_ind = i*2;
        H_mat0_1 = H_mat(a(1), alpha(1), d(1), q(1, start_ind));
        H_mat4_5 = H_mat(a(5), alpha(5), d(5), q(5, start_ind));
        H_mat5_6 = H_mat(a(6), alpha(6), d(6), q(6, start_ind));

        H_mat1_4 = H_mat0_1\robot_end_mat/(H_mat4_5*H_mat5_6);
        P1_3 = H_mat1_4*[0 -link(4) 0 1]' -[0 0 0 1]';
        value = acos((P1_3(1,1)^2 + P1_3(2,1)^2 -link(2)^2 - link(3)^2)/(2*link(2)*link(3)));
        theta3 = real([value -value]);
        
        q(3,start_ind) = mod(theta3(1) + pi, 2*pi) - pi;
        q(3,end_ind) = mod(theta3(2) + pi, 2*pi) - pi;
    %%%%%%%%%%%%%%%%%%%%%% theta 2 %%%%%%%%%%%%%%%%%%%%%%%%
        P1_3_x = P1_3(1,1);
        P1_3_y = P1_3(2,1);

        theta2_1 = -atan2(P1_3_y, -P1_3_x) + asin(link(3)*sin(q(3,start_ind))/(sqrt(P1_3_x^2 + P1_3_y^2)));
        theta2_2 = -atan2(P1_3_y, -P1_3_x) + asin(link(3)*sin(q(3,end_ind))/(sqrt(P1_3_x^2 + P1_3_y^2)));
        theta2 = real([theta2_1 theta2_2]);
        q(2,start_ind) = mod(theta2(1) + pi, 2*pi) - pi;
        q(2,end_ind) = mod(theta2(2) + pi, 2*pi) - pi;
    %%%%%%%%%%%%%%%%%%%%%% theta 4 %%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1 : 2
            theta_start_ind = 2*(i-1) + (j-1) + 1;
            H_mat1_2 = H_mat(a(2), alpha(2), d(2), q(2, theta_start_ind));
            H_mat2_3 = H_mat(a(3), alpha(3), d(3), q(3, theta_start_ind));
            H_mat3_4 = (H_mat1_2*H_mat2_3)\H_mat1_4;
           
            theta4 = real(atan2(H_mat3_4(2,1), H_mat3_4(1,1)));
            q(4, 2*(i-1) + (j-1) + 1) =  mod(theta4 + pi, 2*pi) - pi;
        end
    end
end

function link_mat = H_mat(a, alpha, d , theta)
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    link_mat = [cos_theta, -sin_theta*cos_alpha,   sin_theta*sin_alpha,  a*cos_theta;
                sin_theta,  cos_theta*cos_alpha, -cos_theta*sin_alpha, a*sin_theta;
                0,                 sin_alpha,            cos_alpha,           d;
                0,                      0,                    0,                1];
end

