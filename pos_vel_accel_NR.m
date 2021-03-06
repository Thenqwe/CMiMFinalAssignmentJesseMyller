function [T, Q, QP, QPP] = pos_vel_accel_NR(C_fun, Cq_fun, Ct_fun, Ctt_fun, tend, q_0, dt)

N_t = floor(round(tend/dt));
T = linspace(0, N_t*dt, N_t+1)';
Q = zeros(N_t+1, length(q_0));
QP = zeros(N_t+1, length(q_0));
QPP = zeros(N_t+1, length(q_0));

% [x, iteration_counter] = NR_method(F, J, x, eps)

[qi, icnt] = NR_method(@(q) C_fun(0, q), @(q) Cq_fun(0, q), q_0, 1e-8);
Cqi = Cq_fun(0, qi);
qip = -Cqi\Ct_fun(0, qi);
gr = Ctt_fun(0, qi, qip); % Calculating for acceleration, this may or may
% not be correct as I need it for dynamic analysis aswell, but in this format
% it is not usable there
qipp = Cqi\gr;
Q(1, :) = qi';
QP(1, :) = qip';
QPP(1, :) = qipp';
% Step equations forward in time
for n = 1 : N_t
    [qi, icnt] = NR_method(@(q) C_fun(T(n + 1), q), ...
        @(q) Cq_fun(T(n + 1), q), ...
        qi, 1e-8);
    Cqi = Cq_fun(T(n + 1), qi);
    qip = -Cqi\Ct_fun(T(n + 1), qi);
    gr = Ctt_fun(T(n + 1), qi, qip);
    qipp = Cqi\gr;
    disp(icnt)
    Q(n + 1, :) = qi';
    QP(n + 1, :) = qip';
    QPP(n + 1, :) = qipp';
end