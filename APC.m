%----------------------------------------------------------------------------------
%-------------------         >>>  INITIALIZATIONS  <<<        ---------------------
%----------------------------------------------------------------------------------
ki=0;
NetDef = ['HHHHH';'L----'] ;
NN=[2 2 1];
% ----------   Initializations   -----------
Ts = 0.2;                    % Sampling period (in seconds)
samples = 300;              % Number of samples in simulation
% ------------ Reference filter ---------------
% Am = [1];                    % Denominator of filter
% Bm = [1];                    % Numerator of filter (starts in q^{-1})
apcinit
% ---------- APC initializations -----------
N1 = 1;                      % Minimum prediction horizon (typically=nk)
N2 = 7;                      % Maximum prediction horizon (>= nb)
% Nu = 2;                      % Control horizon
rho = 0.1;                  % Weight factor on differenced control signal

% >>>>>>>>>>>>>>>>>>>>>>>>   DETERMINE REGRESSOR STRUCTURE   <<<<<<<<<<<<<<<<<<<<<<   
na      = NN(1);                       % # of past y's to be used in TDL
nb      = NN(2);                       % # of past u's to be used in TDL
nk      = NN(3);                       % Time delay in system
nab     = na+sum(nb);                  % Number of inputs to each net
outputs = 1;                           % # of outputs is 1 (SISO system)
inputs  = nab;                         % # of inputs
phi     = zeros(inputs,1);             % Initialize regressor vector


% >>>>>>>>>>>>>>>>>    DETERMINE STRUCTURE OF NETWORK MODEL     <<<<<<<<<<<<<<<<<<<
hidden   = length(NetDef(1,:));        % Number of hidden neurons
L_hidden = find(NetDef(1,:)=='L')';    % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';    % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';    % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';    % Location of tanh output neurons
y1       = [zeros(hidden,1)];          % Hidden layer outputs
yhat     = zeros(outputs,1);           % Network output
% y        = 0;

% 1
% 2
% 3


%>>>>>>>>>>>>>>>>>>>>>>>        INITIALIZE VARIABLES        <<<<<<<<<<<<<<<<<<<<<<
% Determine length of reference filter polynomials
% nam = length(Am);
% nbm = length(Bm);

% Initialization of past signals
maxlen = 5;                           
% Initialization of APC polynomials and matrices
% F = zeros(1,na+1);
% E = [1 zeros(1,N2-1)];
% G = zeros(1,N2+nb-1);
% GAMMA   = zeros(N2-nk+1,Nu);
% PHI   = zeros(N2,1);
% F0    = zeros(N2,1);
% rhoI  = rho*eye(Nu);
%invHG = zeros(Nu,N2-N1+1);
% A     = [1 zeros(1,na)];
% B     = zeros(1,nb);
% Miscellaneous initializations
yhat = 0;
% W2=W2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% persistent ref
y_old = zeros(maxlen,1); 
u_old = 0.6*ones(5,1);
difu = 0; 
invMG = zeros(N2,N2-N1+1); 
H = zeros(N2,1); 
S = zeros(N2,3); 
u2 = 0; 

% if isempty(ref), ref=zeros(samples,1);end

% % % %  ref = zeros(samples,1);
% % % % % ref=[7*ones(1000,1);7*ones(1000,1);7*ones(1000,1);...
% % % % %        7*ones(1000,1);7*ones(1000,1)];
% % % % %>>>>>>>>>>>>>>>>    CALCULATE REFERENCE SIGNAL & FILTER IT     <<<<<<<<<<<<<<<<<<
% % % % for jj = 1:samples,
% % % %     ref(jj) = 4*sin(2*pi*50*Ts*(jj-1)-0.5*pi) + 5;
% % % % %       ref(i) = siggener(Ts*(i-1),sq_amp,sq_freq,sin_amp,sin_freq,dc,sqrt(Nvar));     
% % % % end
% fmax=50;
if strcmp(refty,'siggener'),
  ref = zeros(samples+N2,1);
  for i = 1:samples+N2,
    ref(i) = siggener(Ts*(i-1),sq_amp,sq_freq,sin_amp,sin_freq,dc,sqrt(Nvar));
  end
elseif strcmp(refty,'none'),
  ref = zeros(samples+N2,1);
else
  eval(['ref = ' refty ';']);
  ref=ref(:);
  i=length(ref);
  if i>samples+N2,
    ref=ref(1:samples+N2);
  else
    ref=[ref;ref(i)*ones(samples+N2-i,1)];
  end
end
ref=filter(Bm,Am,ref);

t    = -Ts;
yhat = y_0;
%------------------------------------------------------------------------------
%-------------------         >>>   MAIN LOOP   <<<           ------------------
%------------------------------------------------------------------------------
for i=1:samples,
%>>>>>>>>>>>>>>>>>>>>>>>>    READ OUTPUT FROM PLANT     <<<<<<<<<<<<<<<<<<<<<<<
  if strcmp(simul,'simulink')
    utmp=[t-Ts,u_old(1);t,u_old(1)];
    simoptions.InitialState=x0;
    [time,x0,y] = sim(sim_model,[t-Ts t],simoptions,utmp);
    x0 = x0(size(x0,1),:)';
    y  = y(size(y,1),:)';
  elseif strcmp(simul,'matlab')
    ugl = u_old(nk);
    [time,x] = ode45(mat_model,[t-Ts t],x0);
    x0 = x(length(time),:)';
    eval(['y = ' model_out '(x0);']);
  elseif strcmp(simul,'nnet')
    y=yhat;
  end
  
%>>>>>>>>>>>>>>>>>>>>>>     CALCULATE CONTROL SIGNAL     <<<<<<<<<<<<<<<<<<<<<<
e = ref(i) - y;
%>>>>>>>>>>>>>>>>>>>>>>     CALCULATE CONTROL SIGNAL     <<<<<<<<<<<<<<<<<<<<<<
% %  PHI = PHI+F0*y;         % Complete calculation of the vector PHI
U=invMG*(ref(i:i+6)-S*[y;y_old(1);y_old(2)]-H*difu);%invHG*(ref(i+N1-1:i+N2-nk)-PHI(N1:N2));

% if abs(e)<=0.05,
      ucom = u2 + ki*e;
%   else
%       ucom = 0;
% end

u = U(1)+u_old(1)+ucom;
difu = U(1)+ucom;  
 if u>12,
    u=12;
 elseif u<-12
    u=-12;
 end
 %>>>>>>>>>>>>>>>>>>>>>>>>>>       TIME UPDATES        <<<<<<<<<<<<<<<<<<<<<<<<<
  ucon=u;  
  y_old    = shift(y_old,y);
  u_old    = shift(u_old,u);
  u2       = ucom;
%------------------------------------------------------------------------------
%-----------      >>>   DESIGN CONTROLLER FOR NEXT SAMPLE   <<<      ----------
%------------------------------------------------------------------------------
%>>>>>>>>>>>>>>>>>>>>>  CALCULATE OUTPUT PREDICTED BY NN   <<<<<<<<<<<<<<<<<<<<
   phi      = [y_old(1:na);u_old(nk:nk+nb-1)];
   h1 = W1(:,1:inputs)*phi + W1(:,inputs+1);  
   y1(H_hidden) = pmntanh(h1(H_hidden)); 
   y1(L_hidden) = h1(L_hidden);
   h2 = W2(:,1:hidden)*y1 + W2(:,hidden+1);
   yhat(H_output) = pmntanh(h2(H_output));
   yhat(L_output) = h2(L_output);
   if strcmp(simul,'nnet')
     y = yhat;
   end



%>>>>>>>>>>>>>>>>>>>>>>   GET LINEAR PARAMETERS FROM NN  <<<<<<<<<<<<<<<<<<<<<<
   % Matrix consisting of the partial derivatives of each output with
   % respect to each of the outputs from the hidden neurons
   d21 = W2;
%    for j = H_output',
%      d21(j,:) = (1-yhat(j)*yhat(j))*W2(j,:);
%    end

   % Matrix with partial derivatives of the output from each hidden neurons
   % with respect to each input:
   d10 = W1;
   for j = 1:5%H_hidden',
     d10(j,:) = (1-y1(j)*y1(j))*W1(j,:);
   end

   % Matrix with partial derivative of each output with respect to each input
   d20 = d21(1:hidden)*d10;

   A = [1 -d20(1,1:na)];
   Aforlin =-d20(1,1:na);
   B = d20(1,na+1:nab);
   
% ylinear= [-Aforlin,B]*phi; 

%>>>>>>>>>>>>>>>>>>> construct predictor's matrix  <<<<<<<<<<<<<<<<<<<
a1 = Aforlin(1);
a2 = Aforlin(2);
b0 = B(1);
b1 = B(2);

% direct predictor
Gpre = diag(b0*ones(7,1)) + diag(b1*ones(6,1),-1);
Hpre = [b1;0;0;0;0;0;0];
Spre = [1-a1 a1-a2 a2;a1-a2 a2 0;a2 0 0;zeros(4,3)];
Kpre = diag(ones(7,1)) + diag((a1-1)*ones(6,1),-1) + diag((a2-a1)*ones(5,1),-2) + diag(-a2*ones(4,1),-3);
% K = inv(Kpre);
% standard predictor
G = Kpre\Gpre;
H = Kpre\Hpre;
S = Kpre\Spre;
% matrix for control
M = G'*G+rho*eye(7);
invMG = inv(M)*G';
% P0 = -G'*S;
% P1 = -G'*H;
% P2 = -G';
end


% % % %>>>>>>>>>>>>>>>>>>>       SOLVE DIOPHANTINE EQUATIONS     <<<<<<<<<<<<<<<<<<<
% % %   Atilde  = [A 0];                         %              -1
% % %   Atilde(2:na+2)=Atilde(2:na+2)-A(1:na+1); % Atilde = (1-q  )A
% % %   F       = -Atilde(2:na+2);               % F(k) = q(1-Atilde)
% % %   G = [B zeros(1,N2-1)];                   % G(1) = B
% % %   GAMMA(1,1) = G(1);
% % %   F0(1) = F(1);                            % Store F(1)_0
% % %   PHI(1)= G(3-nk:nb)*(u_old(1:nb+nk-2)-u_old(2:nb+nk-1))+...
% % %              F(2:na+1)*y_old(1:na);        % F(k)y(t)-F(k)_0*y(t) 
% % %   for k=1:N2-1,
% % %     G(k+1:k+nb) = G(k+1:k+nb)+B*F(1);      %  G(k+1)    = G(k)+ q^{-1}*B*F(k)_0 
% % %     F      = [F(2:na+1) 0]-F(1)*Atilde(2:na+2); % F(k+1)= q(F(k)-Atilde*F(k)_0)
% % %     PHI(k+1) = G(k-nk+3:k+nb)*(u_old(1:nb+nk-2)-u_old(2:nb+nk-1))+...
% % %              F(2:na+1)*y_old(1:na);        % (G-??)u(t)+F(k)y(t)-F(k)_0*y(t)
% % %     F0(k+1)  = F(1);                       % Store F(k)_0
% % %   end
% % %   
% % %   % Insert G coefficients in GAMMA
% % %   for j=1:Nu,
% % %     GAMMA(j:N2-nk+1,j) = G(1:N2-nk+2-j)';
% % %   end
% % %   H = GAMMA(N1:N2,:)'*GAMMA(N1:N2,:)+rhoI;
% % %   invHG=inv(H)*GAMMA(N1:N2,:)';
 
% end
%------------------------------------------------------------------------------
%------------------        >>>   END OF MAIN LOOP   <<<       -----------------
%------------------------------------------------------------------------------
