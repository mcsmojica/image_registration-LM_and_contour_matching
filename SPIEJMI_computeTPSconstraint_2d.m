function h = SPIEJMI_computeTPSconstraint_2d(v,Tpoints)
    % computes the LHS h(v) of the system of constraints at the current 
    % TPS parameters:
    % n = number of landmarks
    % c_1^i + c_2^i + ... + c_n^i = 0
    % c_1^i*t_1^1 + c_2^i*t_2^1 + ... + c_n^i*t_n^1 = 0
    % c_1^i*t_1^2 + c_2^i*t_2^2 + ... + c_n^i*t_n^2 = 0
    % for i = 1,2
    
    % Inputs:
    % - v = 2*(n+3)x1 vector
    % - Tpoints = nx2 array
    num_LM = size(Tpoints,1);
    v1 = v(1:num_LM); v2 = v(num_LM+4:end-3);
    h = zeros(6,1);
    h(1) = sum(v1);
    h(2) = dot(v1,Tpoints(:,1));
    h(3) = dot(v1,Tpoints(:,2));
    h(4) = sum(v2);
    h(5) = dot(v2,Tpoints(:,1));
    h(6) = dot(v2,Tpoints(:,2));
end