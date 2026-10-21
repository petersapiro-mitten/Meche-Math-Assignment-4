ForwardEuler = struct();
ForwardEuler.A = 0;
ForwardEuler.B = 1;
ForwardEuler.C = 0;


Midpoint = struct();
Midpoint.A = [0 0; 0.5 0];
Midpoint.B = [0 1];
Midpoint.C = [0 0.5];

test_rate = @(t,X) -X;

[XB, num_evals] = explicit_RK_step(test_rate, 0, 1, 0.1, ForwardEuler)


[t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(test_rate, [0 1], 1, 0.1, Midpoint)

X_list(end)
num_evals
exp(-1)