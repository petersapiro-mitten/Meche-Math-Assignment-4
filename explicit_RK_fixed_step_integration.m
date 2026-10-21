%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct)

%find last entry of the tspan
t_end = tspan(end);

%find first entry of the t_span
t_start = tspan(1);

%estimte N
rough_N = (t_end-t_start)/h_ref;

%round N to evenly divide
N = ceil(rough_N);

%find h
h = (t_end-t_start)/N;

%find m rows 
m = length(X0);

%intialize t_list and set first row
t_list = zeros(N + 1,1);
t_list(1) = t_start;

%intialize X_list and set first row
X_list = zeros(N + 1, m);
X_list(1,:) = X0';

%start num evals count
num_evals = 0;

%find h_avg which is just h in this case
h_avg = h;

for k  = 1 : N
    % find current time
    current_t = t_list(k);

    %find current state
    current_state = X_list(k,:)';

    %pass into explicit_Rk_step 
    [XB,step_num_evals] = explicit_RK_step(rate_func_in,current_t,current_state,h,BT_struct);
    
    % store iteration of XB
    X_list(k + 1,:) = XB';

    % store iteration of t
    t_list(k + 1) = current_t + h;

    %update num_evals
    num_evals = num_evals + step_num_evals;
end 
end