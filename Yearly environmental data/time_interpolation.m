function [ x ] = time_interpolation( t_0, x_0 , t_1 , x_1 , t)
% This function evaluates the value of a variable x for a certain date
% based on two values of this variable for dates framing the given one.

%% INPUTS

% t_0: date of the initial value available (d)
% x_0: value of the variable x available for the initial date (any unit)
% t_1: date of the final value available (d)
% x_1: value of the variable x available for the final date (any unit)
% t: date for which the value of x is searched for (d)

%% OUTPUTS

% x: value of the variable computed for the date t (any unit)

%% Calculations

x = x_0 + (x_1 - x_0)/(t_1 - t_0)*(t - t_0);

end

