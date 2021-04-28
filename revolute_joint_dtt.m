function gr = revolute_joint_dtt(i, j, s_i, s_j, q, q2)

idx_i = body_idx(i);
r_i = q(idx_i(1:2));
phi_i = q(idx_i(3));
idx_j = body_idx(j);
r_j = q(idx_j(1:2));
phi_j = q(idx_j(3));

idx_i2 = body_idx(i);
r_i2 = q2(idx_i2(1:2));
phi_i2 = q2(idx_i2(3));
idx_j2 = body_idx(j);
r_j2 = q2(idx_j2(1:2));
phi_j2 = q2(idx_j2(3));

% Calculating gr so acceleration can be calculated for all the bodies 
gr = rot(phi_i) * s_i * phi_i2^2 - rot(phi_j) * s_j * phi_j2^2;