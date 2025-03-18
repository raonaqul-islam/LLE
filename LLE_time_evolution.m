%% LLE time evolution
% ----------------------------+
% Raonaqul Islam, UMBC        |
% Date started: March 2, 2025 |
% Last updated: March 9, 2025 |
% ----------------------------+

clear
clc
close

% Parameters
% -----------
alpha   = 3.5;                            % Detuning
beta    = -0.1549;                        % Dispersion, beta2 only in this case
F       = 2.3;                            % Power
gamma   = 1;                              % Normalized nonlinear coefficient

% Spatial domain discretization
% ------------------------------
N       = 512;                            % No. of points on axis
naxis   = (-N/2:N/2-1).';                 % General axis
dtheta  = 2*pi/N;                         % Spatial step-size
theta   = naxis*dtheta;                   % Spatial domain
dmu     = 1;                              % Mode number domain step-size
mu      = fftshift(dmu*naxis);            % Mode number domain

% Power and dispersion tuning
% ----------------------------
F_tilde = fft(F*ones(1,N)).';             % Fourier transform of power

% Time domain discretization
% ---------------------------
dt      = 1e-4;                           % Temporal step size
N_tau_p = 100;                            % No. of photon lifetimes
N_steps = N_tau_p/dt;                     % No. of steps

% Generate figure to plot on
figure('Units','inches','Position',[6 4 8 4]); 

% Time evolution
% ---------------
psi_0   = 2*sech(5*theta);                % Initial guess

for ix = 1:N_steps

    % Nonlinear part (half-step)        
    psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

    % Linear part (full-step)
    psi_0_tilde = fft(psi_nl);    
    A_tilde     = -(1+1i*alpha)+1i.*beta.*(mu).^2;      
    psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
    psi_l       = ifft(psi_l_tilde);
    
    % Nonlinear part (half-step)
    psi_0       = psi_l;
    psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
    psi_0       = psi_nl2; 

    % Output after one-step
    psi_out     = psi_nl2;
    
    % Plot solution for every photon lifetime (1/dt)
    if rem(ix,1/dt) == 0
        plot(theta,abs(psi_out).^2,'LineWidth',2.5)
        xlabel('\theta')
        ylabel('|\psi|^2')
        title(['No. of \tau_{p} = ' num2str(ix*dt)])
        set(gca,'FontSize',16)  
        xticklabels(strrep(xticklabels,'-','−'))
        grid on
        drawnow
    end

end

% Plot soliton in mode number domain
% -----------------------------------
% Prepare the solution
psi_fft  = abs(fftshift(fft(psi_out))).^2;  % Obtain fourier transform of the solution
psi_norm = psi_fft/max(psi_fft);            % Normalize the output
psi_db   = 10*log10(psi_norm);              % Calculate normalized output power in dB

% Plot in stem
figure('Units','inches','Position',[6 4 15 4]);
soliton = stem(fftshift(mu),psi_db,'LineWidth',2,'Marker','None','Color','Red');
soliton.BaseValue = -400; % Turn the plot upside down
xlabel('Relative mode number, \mu')
ylabel('Power (dB)')
set(gca,'FontSize',16)  
xticklabels(strrep(xticklabels,'-','−'))
yticklabels(strrep(yticklabels,'-','−'))
grid on