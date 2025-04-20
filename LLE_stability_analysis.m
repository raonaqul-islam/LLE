% LLE stability analysis with fsolve
% -----------------------------+
% Raonaqul Islam, UMBC         |
% Date started: March 22, 2025 |
% Last updated: April 20, 2025 |
% -----------------------------+

clear
close

% Input Parameters
% -----------------
alpha   = 3.5;                                  % Detuning
beta    = -0.1549;                              % Dispersion, beta2 only in this case
F       = 2.3;                                  % Power
gamma   = 1;                                    % Normalized nonlinear coefficient
method  = 'fft';                                % Method for differentiation, 'fdm' or 'fft'

% Spatial domain discretization
% ------------------------------
N       = 512;                                  % No. of points on axis
naxis   = (-N/2:N/2-1).';                       % General axis
dtheta  = 2*pi/N;                               % Spatial step-size
theta   = naxis*dtheta;                         % Spatial domain
dmu     = 1;                                    % Mode number domain step-size
mu      = fftshift(dmu*naxis);                  % Mode number domain

% fsolve operation
% -----------------
% Initial guess
psi_r0  = 2*sech(5*theta);                      % A common guess for solitons
psi_m0  = zeros(N,1);                           % Initially, no imaginary parts in the guess
psi_0   = [psi_r0;psi_m0];                      % Complete guess

% fsolve function
switch method
    case 'fdm'
        fun     = LLE_fdm(alpha,beta,gamma,F,dtheta,mu,N); 
    case 'fft'
        fun     = LLE_fft(alpha,beta,gamma,F,dtheta,mu,N);
end

% fsolve options
options = optimset('Display','iter','Jacobian','on','TolFun',1e-12,'TolX',...
          1e-12,'Algorithm','levenberg-marquardt','ScaleProblem','Jacobian',...
          'MaxIter',100);

% find roots using fsolve
[psi,fval,exitflag,output,Jacobian] = fsolve(@fun.findroots,psi_0,options);    

% Extract and form the complex solution
% --------------------------------------
psi_r   = psi(1:N);
psi_m   = psi(N+1:end);
psi_out = psi_r + 1i*psi_m;

% Plot Soliton
% -------------
figure('Units','Inches','Position',[2 2 12 8]);
subplot(211)
plt1    = plot(theta./pi,abs(psi_out),'Color',[0.07 0.62 1 0.8]);
xtext   = '\theta/\pi';
ytext   = '|\Psi|';
customplot(plt1,xtext,ytext);

% Plot Eigenvalues
% -----------------
eigenvalues = eig(Jacobian);
subplot(212)
plt2     = scatter(real(eigenvalues),imag(eigenvalues),'MarkerEdgeColor',[0.07 0.62 1],'MarkerFaceColor',[0.07 0.62 1]);
xtext    = 'Real';
ytext    = 'Imaginary';
xline(0,'LineWidth',3,'LineStyle','--','Color','k');
customplot(plt2,xtext,ytext);
