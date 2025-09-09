    
    %=================================================================
    %=================================================================
    
    % "A variational and symplectic framework for model-free control: preliminary results"
    
    % Code associated to the work submitted to the CCTA'25 conference
    
    % (c) [2025]  Nantes Université - Centrale Nantes - LS2N UMR 6004, Nantes
    % by Loïc MICHEL
    % All rights reserved under MIT license.
    
    % The file 'ModelFreeControl_Symplectic_PwrConverter_SweepEpsilon.m' computes both the evolution of gamma(t) 
    % and the tracking error with respect to several values of Epsilon_M (Figs. 12 and 13).
    
    %=================================================================
    %=================================================================
    
    
    clc
    clear all
    close all
    warning off
    
    %---------------------------------------
    PowerControl = 0;
    % '0' for voltage control
    % '1' for power control
    VarIntegrator = 1; %should not be changed
    % '0' for standard model-free control
    % '1' for variationnal-based model-free control
    %---------------------------------------
    
    
  % Simulation parameters
    tspan = [0 0.01];
    N = 10000;
    h = 1e-6;
  
    % Buck circuit - nominal parameters
    
    R = 10;
    L = 1e-3;
    C = 10e-6;
    Vin = 20;
    
    % Disturbance management
    
    FirstPerturbation = 1;
    TimeFirstPerturbation = 0.004;
    
    SecondPerturbation = 1;
    TimeSecPerturbation = 0.006;
    
    Delta_R = 5;   % increase of resistance - unstable perturbation (voltage cont.)
    %Delta_R = 1;   % increase of resistance - stable perturbation (voltage cont.)
    Delta_C = 1e-9;   % increase of capacitor

    saturation = 1;            % Enable saturation on the control output 

         
    if PowerControl == 0

        % ------- Voltage control settings -------
    
        % Tracking ref. settings
        const = 2.5;
        type_consigne = 2;
        slope = 1000;
        t_start_slope = 0.005;
    
        alpha = 10; % model-free control key parameter
        gamma0 = 1; % init. of the variational integrator
        Kp = 2; % proportional corrector (in the model-free law)
    
        C_gamma = 1; % variational param
        Gamma_ = h^2 / C_gamma;
        Gamma_0 = Gamma_;
        Epsilon_M = 0.1; % variational param
    
        % ----------------------------------------

    else
    
        % ------- Power control settings -------

        % Tracking ref. settings
        const = 0.25;
        type_consigne = 1;
        slope = 100;
        t_start_slope = 0.005;
    
        %Control parameters
        alpha = 3; % model-free control key parameter
        gamma0 = 3; % init. of the variational integrator
        Kp = 2; % proportional corrector (in the model-free law)

        C_gamma = 1; % variational param
        Gamma_ = h^2 / C_gamma;
        Epsilon_M = 0.1; % variational param

        % ----------------------------------------
    end
        
    t = tspan(1)+[0:N]*(tspan(2) - tspan(1))/N;
    
    
    hh = 1e-3; % factor for the computation of the time-derivatives ->
    % must be integrated within the alpha and Gamma parameters for more
    % visibility ;)

    EpsilonM_vec = [Epsilon_M/1e4, Epsilon_M, Epsilon_M*1e4 ];

    for IndexEpsilon = 1:3

        Epsilon_ = EpsilonM_vec(IndexEpsilon);

        fprintf('Epsilon_M = %e \n', Epsilon_);
    
        % Simulation loop
        for k = 1:N
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% Compute the reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        if PowerControl == 1

        % A first slope, then a constant
            if slope*t(k) <= const
                y_reference(k) = slope*t(k);
            else
                y_reference(k) = const;
            end
    
        else 

        % A first slope, then a constant and afterwards a slope
            if slope*t(k) <= const
                y_reference(k) = slope*t(k);
            else
                y_reference(k) = const;
            end
    
            if t(k) >= t_start_slope
                y_reference(k) = const + slope*(t(k) - t_start_slope)  ;
            end
    
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Computation of the first order derivatives %%%%%%%%%%%%%%%%%%%%
 
        % Computation of the first order derivative for the reference
        if k > 1
            dy_reference(k-1) = (y_reference(k) - y_reference(k-1))/hh;
        else
            dy_reference(k) = 0;
        end
    
        % Computation of the first order derivative for the output y   
        if k > 1
            dy(k-1) = (y(k) - y(k-1))/hh;
        else
            dy(k) = 0;
        end
    
        % Computation of the P-action
        if k > 1
            P_action(k) = Kp * (y_reference(k) - y(k-1));
        else
            P_action(k) = 1;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Computation of the variational integrator + model-free control -------
        if VarIntegrator == 1
    
            if k > 1
    
                %%%%
                Depsilon = (dy(k-1) - dy_reference(k-1) );
        
                if k > 3
        
                    % Computation of the variational integrator
                    % (Epsilon_ is modified w.r.t. the 'EpsilonM_vec' vector - Gamma_0 is kept constant)
                    gamma_t(k) = (Depsilon - Epsilon_)*Gamma_0 + 2*gamma_t(k-1) - gamma_t(k-1);
        
                    gamma_t(k-2) = gamma_t(k-1);
                    gamma_t(k-1) = gamma_t(k);
    
                else
    
                    gamma_t(k) = gamma0;
    
                end
                
                % Computation of the variationnal-based model-free control law
                u_(k) = u_(k-1) - (gamma_t(k))*(dy(k-1) - dy_reference(k-1)) +  P_action(k) - C_gamma * ((gamma_t(k) - gamma_t(k-1))/hh); 
    
            else

                u_(k) = P_action(k);
                gamma_t = gamma0;
    
            end
    
        end
        % --------------------------------------------------------------  

        % ------- Computation of the standard model-free control -------
        if VarIntegrator == 0
    
            if k > 1
                u_(k) = (u_(k-1) -1/alpha*(dy(k-1) - dy_reference(k-1)) +  P_action(k))/1;  %
            else
                u_(k) = P_action(k);
            end
    
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%% Ouput of the control management %%%%%%%%%%%%%%%%%%%%%%%

        outputGain = 2000;   % output gain
        duty_cycle = 1/outputGain * u_(k);
    
    
        if saturation == 1
    
            if duty_cycle > 0.9
                duty_cycle = 0.9;
            elseif duty_cycle < 1e-3
                duty_cycle = 1e-3;
            end
    
        end
    
        duty_cycle_vec(k) = duty_cycle;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
        % Perturbation management (needs 'FirstPerturbation' or
        % 'SecondPerturbation' enabled)
        if t(k) >= TimeFirstPerturbation && FirstPerturbation == 1
            FirstPerturbation = 0;
            R = R + Delta_R;
            C = C + Delta_C;
        end
    
        if t(k) >= TimeSecPerturbation && SecondPerturbation == 1
            SecondPerturbation = 0;
            R = R + Delta_R;
            C = C + Delta_C;
        end
        
        % Update of the buck averaged model
        A = [0 -1/L;1/C -1/(R*C)];
        u = [duty_cycle/L*Vin; 0];
        B = [1  0];
        
        %Runge-Kutta method
        if k == 1
            yRK = zeros(1,length(A))';
            y(1) = 0;
        end
    
        f1 = h*av_SS(A, B, u, t(k),yRK(:,k));
        f2 = h*av_SS(A, B, u, t(k) + h/2,yRK(:,k) + f1/2);
        f3 = h*av_SS(A, B, u, t(k) + h/2,yRK(:,k) + f2/2);
        f4 = h*av_SS(A, B, u, t(k) + h,yRK(:,k) + f3);
        yRK(:,k + 1) = yRK(:,k) + (f1 + 2*(f2 + f3) + f4)/6;
    
    
        % Select the model output ('PowerControl' and 'VarIntegrator')
        if PowerControl == 0
            y(k + 1) = (yRK(2,k + 1));
        else
            y(k + 1) = (yRK(2,k + 1)) * (yRK(1,k + 1)) ;
        end
        
      
        end

        tracking_error(IndexEpsilon,:) = y(1:end-1) - y_reference(1:end);

        gamma_curves(IndexEpsilon,:) = gamma_t;

    end


    
    FtSize = 50;
    
    %%%%%%%%% PLOT SECTION
    
    color_vec = ['b', '--r', 'oc', '.g', 'y'];
    
            % figure('name','Gamma variation vs. epsilon_M')
            % for ii = 1:5
            % semilogy(t(1:end-1),  abs(gamma_curves(ii,:)), color_vec(ii), 'linewidth', 3)
            % hold on
            % end
            % grid on
            % xlabel('time (s)','FontSize', FtSize, 'Interpreter','latex');
            % ylabel('$\gamma(t)$','FontSize',FtSize, 'Interpreter','latex');
            % legend('$\varepsilon_M/100$','$\varepsilon_M/10$','$\varepsilon_M$','$10 \varepsilon_M$','$100 \varepsilon_M$', 'FontSize', FtSize, 'Interpreter','latex')
            % set(gcf,'color',[1 1 1]);
            % set(gca,'fontsize',FtSize);
            % ax = gca;
            % ax.YAxis.Exponent = -5;
    

            figure('name','Tracking error vs. epsilon_M')
            for ii = 1:3
            semilogy(t(1:end-1),  abs(tracking_error(ii,:)), color_vec(ii), 'linewidth', 4)
            hold on
            end
            grid on
            xlabel('Time (s)','FontSize', FtSize, 'Interpreter','latex');
            ylabel('tracking error','FontSize', FtSize, 'Interpreter','latex');
            legend('$\varepsilon_M/10^4$','$\varepsilon_M$','$10^4 \varepsilon_M$', 'FontSize', FtSize, 'Interpreter','latex')
            set(gcf,'color',[1 1 1]);
            set(gca,'fontsize', FtSize);
            ax = gca;
            ax.YAxis.Exponent = -5;
    
 

    fprintf("End of the program \n");