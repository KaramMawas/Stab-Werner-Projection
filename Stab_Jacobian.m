function [J] = Stab_Jacobian(LAMDA,PHI,R)

Lamda = LAMDA*pi/180 ;
Phi = PHI*pi/180 ;
delta = (pi/2)-Phi ;
t = Lamda .*(cos(Phi)./delta) ;


J = [-R*sin(t).*cos(Phi) -R*cos(t)+R*Lamda.*sin(t).*sin(Phi)-(R*Lamda./delta).*sin(t).*cos(Phi);
    R*cos(Phi).*cos(t) -R*sin(t)-R*Lamda.*cos(t).*cos(Phi)+(R*Lamda./delta).*cos(t).*cos(Phi)];
