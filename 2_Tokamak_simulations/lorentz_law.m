%% función que calcula la aceleración de la particula

function [dR] = lorentz_law(t,R)

global q m

%% Se extraen r y v de R

r = R(1:3);
v = R(4:6);

%% Se evalua el valor de B y E

B = eval_B (r) ; 
% E = eval_E (r) ;
E=0;

%% Se calcula la aceleración

a = q/m*(E+cross(v,B));

%% Se construye el vector dR

dR = [v;a];
 

% t
end