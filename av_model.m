
%=====================================
%     modèle moyen généralisé
%
%           4 oct 2009
%
%=====================================

clear all
close all
clc

% [0] J. Mahdavi, "Analysis of power electronic converters using the generalized
% state-space averaging approach", 1997.

% paramètres circuit

R = 10;
L = 1e-3;
C = 10e-6;
Ts = 0.1e-3;
w = 2*pi*1/Ts;
d = 0.25;
Vin = 20;

% modèle moyen

% u fermé

% d il / d t = (E - v0)/L

% d v0 / dt = 1/C( il - v0 /R)

% u ouvert

% d il / dt = v0

% d v0 / dt = 1/C ( v0 / R)

% fusion topologique

% d il / d t = 1/L ( E*u - vo)

% d v0 / dt = 1/C (il(1-u) - v0/R)

% => il(1-u) - I = 0 => 

% => LC d^2 v0 / d t^2 + L/R d v0 / d t = E

% modèle moyen

    % paramètre simulation

    tspan = [0 0.01];
    N = 10000;

A = [0 -1/L;1/C -1/(R*C)];
B = [1  1];

u = [d/L*Vin; 0];

        %Runge-Kutta method to solve vector differential eqn y’(t) = f(t,y(t))
        % for tspan = [t0,tf] and with the initial value y0 and N time steps
         
            h = (tspan(2) - tspan(1))/N; %t = tspan(1)+[0:N]'*h;
       
            t = tspan(1)+[0:N]*(tspan(2) - tspan(1))/N;

            y = zeros(1,length(A))';

                for k = 1:N
            
            	%u = input(t);
            
                    f1 = h*av_SS(A, B, u, t(k),y(:,k));
                    f2 = h*av_SS(A, B, u, t(k) + h/2,y(:,k) + f1/2); 
                    f3 = h*av_SS(A, B, u, t(k) + h/2,y(:,k) + f2/2); 
                    f4 = h*av_SS(A, B, u, t(k) + h,y(:,k) + f3);
                    y(:,k + 1) = y(:,k) + (f1 + 2*(f2 + f3) + f4)/6;
                end

figure('Name','Modèle moyen','NumberTitle','on')
plot(t,y,'r')
legend('courant moyen')

ymoyen = y;


h = (tspan(2) - tspan(1))/N; %t = tspan(1)+[0:N]'*h;
       
            t = tspan(1)+[0:N]*(tspan(2) - tspan(1))/N;

           A = [0 w -1/L 0 0 0;-w 0 0 -1/L 0 0;1/C 0 -1/(R*C) w 0 0; 0 1/C -w -1/(R*C) 0 0; 0 0 0 0 0 -1/L; 0 0 0 0 1/C -1/(R*C)];
           B = [1  1 1  1  1  1];
                            
            yg = zeros(1,length(A))';

               for k = 1:N

                         d = 0.25; %u(k);  

                            if d >= 0.9
                            d = 0.9;
                            elseif d <= 0.1
                            d = 0.1;
                            end
                            

        
                            u_ = [sin(2*pi*d)/(2*pi*L)*Vin ; (sin(pi*d))^2/(pi*L)*(-Vin); 0;0; d/L*Vin; 0];
                                                
                    % Évaluatiok du modèle
                       
                    f1 = h*av_SS(A, B, u_, t(k), yg(:,k));
                    f2 = h*av_SS(A, B, u_, t(k) + h/2,yg(:,k) + f1/2); 
                    f3 = h*av_SS(A, B, u_, t(k) + h/2,yg(:,k) + f2/2); 
                    f4 = h*av_SS(A, B, u_, t(k) + h,yg(:,k) + f3);
                    yg(:,k + 1) = yg(:,k) + (f1 + 2*(f2 + f3) + f4)/6;
                    
                   y(k+1) = yg(6,k+1)+ 2.*yg(3,k+1).*cos(w*t(k+1))- 2.*yg(4,k+1).*sin(w*t(k+1));
                   
                end
                               
%PSIMtime = 0;
%PSIMvout = 0;
%[PSIMtime,PSIMvout] = importPSIM('PWM_hacheur.txt',PSIMtime,PSIMvout,3);

%plot(t(1:ekd-1), yC);

figure('name','Modèle moyen généralisé complexe')
subplot(2,1,1)
plot(t, yg(5,:)+ 2.*yg(1,:).*cos(w*t)- 2.*yg(2,:).*sin(w*t))
legend('courant moyen ordre zéro','courant moyen','courant moyen généralisé');
subplot(2,1,2)
plot(t,yg(6,:)+ 2.*yg(3,:).*cos(w*t)- 2.*yg(4,:).*sin(w*t))
legend('teksiok moyenne généralisé', 'simulation PSIM');



figure(3)
plot(t,ymoyen(2,:),'r')
hold on
plot(t,yg(6,:)+ 2.*yg(3,:).*cos(w*t)- 2.*yg(4,:).*sin(w*t))
legend('avaraged model created by Sigma\Phi-PES', 'Equivalent generalized avarage model')
xlabel('time (s)');
ylabel('Vs (V)');





