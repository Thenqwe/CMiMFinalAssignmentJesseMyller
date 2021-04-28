function acc = system_accelerations(t, q, qp, M, sforce, grav, bodies, Cq, C, g)

alfa = 10;
beta = 10;

F = force_vector(grav, sforce, bodies, q);
MCq = [M, transpose(Cq);
    Cq, 0.*Cq];
Fg = [F;
    g - 2 * alfa * 0 - beta^2 * C];
acclambda = MCq\Fg; % This matrix holds values of both accelerations and Lagrange multipliers
acc = acclambda(1:length(acclambda)/2); % We are only interested in the accelerations