function F=duffin(t,y,PAR) 
%%
% This script defines the system to integrate in the ode45 command. In this
% case the system is a duffin oscillator.
% The first two input variables are
% respectively the time and the time variable. The last variable is a
% vector that we use to input the parameters of the oscillator. In the
% variable PAR the components are defined as follows: 
%
% * PAR(1): mass (m)
% * PAR(2): damping (c)
% * PAR(3): linear stiffness (k) 
% * PAR(4): nonlinear stiffness (\alpha)
% * PAR(5): amplitude of excitation
% * PAR(6): frequency of excitation 
%%
% The function has to be written in the form $\dot{x}=F(x,\dot{x})$. The
% output variable F, hence is $\dot{x}$. In this case we use 
%%
% $$y(1)=\dot{x}$$
%
% $$y(2)=x$$
%%
% The ODE $m\ddot{x}+c\dot{x}+kx+\alpha x^3=A\sin(\Omega t)$ becomes 
%%
% $$F(1)=-\frac{c}{m}y(1)-\frac{k}{m}y(2)-\frac{\alpha}{m}y(2)^3+\frac{A}{m}\sin(\Omega t)$$
%
% $$F(2)=y(1)$$

m      =PAR(1);
c      =PAR(2);
k      =PAR(3);
alpha  =PAR(4);% zero
A      =PAR(5);
Omega  =PAR(6);


F(1)=(c/m.*(y(1)^3))+(-c/m.*(y(1)^2))+(-c/m.*y(1))-(k/m.* y(2))-(A*sin(Omega.*t));
F(2)=y(1);
F=F';
end %The function ends here
