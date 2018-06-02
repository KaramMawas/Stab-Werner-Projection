function [x,y]=Stab(LAMDA,PHI,R)
% Name, given name, matrikulation number
% MAWAS,Karam, NUM:2946939
% salih Elankah matrik. nr. 2928326
% calculates the Stab-Werner map projection

Lamda = LAMDA*pi/180 ;
Phi = PHI*pi/180 ;
delta = (pi/2)-Phi ;
t = Lamda .*(cos(Phi)./delta) ;


x=R*delta.*cos(t);
y=R*delta.*sin(t);

