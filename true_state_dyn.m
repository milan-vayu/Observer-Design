% true state dynamics
function xdot=true_state_dyn(t,x,u)
global A B D;
f=[0.6*x(1)*sin(2*t);0.6*x(2)*cos(2*t); 0];
v=2*sin(3*t);
xdot = A*x +B*u + D*v + f;
end