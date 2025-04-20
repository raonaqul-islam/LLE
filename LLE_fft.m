% The only difference this one and the fdm one is that this one uses 
% fourier transform of the operator d^2/dtheta^2.

classdef LLE_fft
    
    % #####################################################################
    properties
        alpha
        beta
        gamma
        F
        dx         
        mu
        N
    end
    
    % #####################################################################
    methods
        
        % #################################################################
        function self = LLE_fft(alpha,beta,gamma,F,dx,mu,N)
            self.alpha = alpha;
            self.beta  = beta;
            self.gamma = gamma;
            self.F     = F;
            self.dx    = dx;            
            self.mu    = mu;
            self.N     = N;
        end

        % #################################################################
        function y = der2(self)            
            y = real(ifft(diag(-self.mu.^2)*fft(eye(self.N))));
        end

        % #################################################################
        function [psi,J] = findroots(self, psi_0)
            psi_r = psi_0(1:self.N,1);
            psi_m = psi_0(self.N+1:end,1);

            F1 = - psi_r + self.alpha.*psi_m + (self.beta/2).*self.der2*psi_m ...
                - self.gamma .* psi_r.^2 .* psi_m - self.gamma .* psi_m.^3 + self.F;

            F2 = - psi_m - self.alpha.*psi_r - (self.beta/2).*self.der2*psi_r ...
                + self.gamma .* psi_r.^3 + self.gamma .* psi_m.^2 .* psi_r;

            psi = [F1;
                   F2];

            J = self.jacobian(psi_0);  % Use psi_0 here
        end

        % #################################################################
        function [J] = jacobian(self,psi_0)
            psi_r = psi_0(1:self.N,1);
            psi_m = psi_0(self.N+1:end,1);

            D2 = self.der2;

            J11 = - eye(self.N) - diag(2 .* psi_r .* psi_m * self.gamma);
            J12 = self.alpha .* eye(self.N) + (self.beta/2)*D2 ...
                - diag(self.gamma .* psi_r.^2) - diag(3 * self.gamma .* psi_m.^2);

            J21 = -self.alpha .* eye(self.N) - (self.beta/2)*D2 ...
                + diag(3 * self.gamma .* psi_r.^2) + diag(self.gamma .* psi_m.^2);

            J22 = - eye(self.N) + diag(2 .* psi_m .* psi_r * self.gamma);

            J = [J11, J12;
                J21, J22];
        end

    end
end