%Field line tracing function

function [dr] = field_line_eq(t,r)

% Se evalua el valor de B

B = eval_B (r) ; 

%% Se construye el vector dR

dr = B;
 

% t
end