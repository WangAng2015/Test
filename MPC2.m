%¶ª·¬Í¼
A = [1 1 2];
B = [3 4];
na =2;
nb =2;
nk =1;
N1 = 1;                      % Minimum prediction horizon (typically=nk)
N2 = 7;                      % Maximum prediction horizon (>= nb)
Nu = 7;                      % Control horizon
rhoI = 0.1;                  % Weight factor on differenced control signal
u_old = [1 2]';
y = 1;
y_old = [3 4]';
PHI   = zeros(N2,1);
F0    = zeros(N2,1);
ref = ones(7,1);
%>>>>>>>>>>>>>>>>>>>       SOLVE DIOPHANTINE EQUATIONS     <<<<<<<<<<<<<<<<<<<
  Atilde  = [A 0];                         %              -1
  Atilde(2:na+2)=Atilde(2:na+2)-A(1:na+1); % Atilde = (1-q  )A
  F       = -Atilde(2:na+2);               % F(k) = q(1-Atilde)
  G = [B zeros(1,N2-1)];                   % G(1) = B
  GAMMA(1,1) = G(1);
  F0(1) = F(1);                            % Store F(1)_0
  PHI(1)= G(3-nk:nb)*(u_old(1:nb+nk-2)-u_old(2:nb+nk-1))+...
             F(2:na+1)*y_old(1:na);        % F(k)y(t)-F(k)_0*y(t) 
  for k=1:N2-1,
    G(k+1:k+nb) = G(k+1:k+nb)+B*F(1);      %  G(k+1)    = G(k)+ q^{-1}*B*F(k)_0 
    F      = [F(2:na+1) 0]-F(1)*Atilde(2:na+2); % F(k+1)= q(F(k)-Atilde*F(k)_0)
    PHI(k+1) = G(k-nk+3:k+nb)*(u_old(1:nb+nk-2)-u_old(2:nb+nk-1))+...
             F(2:na+1)*y_old(1:na);        % (G-??)u(t)+F(k)y(t)-F(k)_0*y(t)
    F0(k+1)  = F(1);                       % Store F(k)_0
  end
  
  % Insert G coefficients in GAMMA
  for j=1:Nu,
    GAMMA(j:N2-nk+1,j) = G(1:N2-nk+2-j)';
  end
  H = GAMMA(N1:N2,:)'*GAMMA(N1:N2,:)+rhoI;
  invHG=inv(H)*GAMMA(N1:N2,:)'
  PHI = PHI+F0*y;         % Complete calculation of the vector PHI
  U=invHG*(ref-PHI(N1:N2))