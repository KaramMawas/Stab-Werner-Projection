close all;
clear all;
clc;

% program settings
R = 6378;
scale = 700;
Phi_intervall    = -90:30:90;
Lambda_intervall = -180:30:180;
Phi_intervall_Tissot    = -60:30:60;
Lambda_intervall_Tissot = -180:30:180;

% open new figure with white background
fig = figure('color', [1 1 1]);

% coastlines
load('coast.mat');
[x,y]= Stab(long,lat,R);
plot(x,y)
hold on

%maridian
[Lamda,Phi]=meshgrid(Lambda_intervall, min(Phi_intervall):1:max(Phi_intervall));
[x,y]=Stab(Lamda,Phi,R);
plot(x,y,'color',[0.5 0.5 0.5]) 
hold on

%parallel circles
[Lamda,Phi]=meshgrid(min(Lambda_intervall):1:max(Lambda_intervall), Phi_intervall);
[x,y]=Stab(Lamda,Phi,R);
plot(x',y','color',[0.5 0.5 0.5])

% Tissot ellipses
for Lamda = Lambda_intervall_Tissot
    for Phi = Phi_intervall_Tissot
        [x, y] = Stab(Lamda, Phi, R);
        
        % metric matrix of the source
        G = R^2 * [cos(Phi * pi / 180)^2 0; 0 1];
        
        % Jacobian
        J = Stab_Jacobian(Lamda, Phi, R);
        
        % Cauchy-Green tensor
        C = J' * J;
        
        % solve the general eigenvalue problem
        [F, Lambda_12] = eig(C, G);
        
        % transform the eigenvectors from source to map
        f = J * F;
        
        % length of the semi axes
        lambda1 = sqrt(Lambda_12(1, 1));
        lambda2 = sqrt(Lambda_12(2, 2));
        
        % angle of semi major axis
        ang = atan2(f(2, 1), f(1, 1));
        
        % ATTENTION!
        % consider the following three cases:
        % 1.) no distortion: green circle of radius r = lambda1 = lambda2 == (equal) 1
        % 2.) conformal:     red circle of radius r = lambda1 = lambda2 ~= (not equal) 1
        % 3.) else:          red ellipse with semi major axis a = lambda1 and semi minor axis b = lambda2
        if abs(lambda1-lambda2)<1e-6
        ellipse(scale * lambda1, scale * lambda2, ang, x, y, 'g');
        elseif  lambda1==lambda2
        ellipse(scale * lambda1, scale * lambda2, ang, x, y, 'r');
        else 
        ellipse(scale * lambda1, scale * lambda2, ang, x, y, 'r');
        end
        
        ca = cos(ang);
        sa = sin(ang);
        
        ax(1, 1) = x - scale * ca * lambda1;
        ay(1, 1) = y - scale * sa * lambda1;
        ax(2, 1) = x + scale * ca * lambda1;
        ay(2, 1) = y + scale * sa * lambda1;

        ax(1, 2) = x + scale * sa * lambda2;
        ay(1, 2) = y - scale * ca * lambda2;
        ax(2, 2) = x - scale * sa * lambda2;
        ay(2, 2) = y + scale * ca * lambda2;
       
        if abs(lambda1-lambda2)<1e-6
        plot(ax, ay, 'g');
        elseif lambda1==lambda2
            plot(ax, ay, 'r');
        else
            plot(ax, ay, 'r');
        end
    end
end

title ('Stab-Werner Projection worded by Karam MAWAS(2946939) Und Salih ElAnkah(2928326)');
axis equal;
axis off;
