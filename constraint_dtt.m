function C = constraint_dtt(revolute, simple, driving, t, q, q2)

r_len = length(revolute);
s_len = length(simple);
d_len = length(driving);

n_constr = 2 * r_len + s_len + d_len;

C = zeros(n_constr, 1);

c_idx = 0;

for r = revolute
    C(c_idx + (1:2)) = revolute_joint_dtt(r.i, r.j, r.s_i, r.s_j, q, q2);
    c_idx = c_idx + 2;
end

for d = driving
    c_idx = c_idx + 1;
    C(c_idx) = driving_joint_dt(d.d_k_t, t); % The result is basically same as with velocity
end