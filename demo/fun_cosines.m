%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : "Cosines" test function
%--------------------------------------------------------------------------
function y = fun_cosines(x1, x2)
  u = 1.6*x1 - 0.5;
  v = 1.6*x2 - 0.5;
  y = 1 - (u.^2 + v.^2 - 0.3*cos(3*pi*u) - 0.3*cos(3*pi*v));
end