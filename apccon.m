% ----------------------------------   APCCON.M    --------------------------------
%
% Program for simulating approximate predictive controllers (APC)
%
% All design parameters must be defined in the file 'apcinit.m'
%
%
% Written by Magnus Norgaard, IAU, Technical University of Denmark.
% LastEditDate: Jan. 23, 2000

%----------------------------------------------------------------------------------
%-------------------         >>>  INITIALIZATIONS  <<<        ---------------------
%----------------------------------------------------------------------------------

%>>>>>>>>>>>>>>>>>>>>>>      READ VARIABLES FROM FILE       <<<<<<<<<<<<<<<<<<<<<<<
clear plot_a plot_b
global ugl
apcinit
eval(['load ' nnfile]);                % Load neural network
NetDef = ['HHHHH';'L----'] ;
NN=[2 2 1];

% >>>>>>>>>>>>>>>>>>>>>>>>   DETERMINE REGRESSOR STRUCTURE   <<<<<<<<<<<<<<<<<<<<<<   
na      = NN(1);                       % # of past y's entering the TDL
nb      = NN(2);                       % # of past u's entering the TDL
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


%>>>>>>>>>>>>>>>>>>>>>>>        INITIALIZE VARIABLES        <<<<<<<<<<<<<<<<<<<<<<
% Determine length of reference filter polynomials
nam = length(Am);
nbm = length(Bm);

% Initialization of past signals
maxlength = 5;                          % MIGHT BE NECESSARY TO INCREASE maxlength
ref_old   = repmat(y_0,maxlength,1);    % FOR HIGH-ORDER SYSTEMS!!!
y_old     = repmat(y_0,maxlength,1);
yhat_old  = repmat(y_0,maxlength,1);
u_old     = repmat(u_0,maxlength,1);

% Initialization of PID parameters
if strcmp(regty,'pid'),
  B1 = K*(1+Ts*Wi/2);
  A1 = Ts*Wi;
  B2 = (2*Td+Ts)/(2*alf*Td+Ts);
  A2 = 2*Ts/(2*alf*Td+Ts);
  I1 = 0;
  I2 = 0;
  uimin = -10; uimax = 10;
end

% Initialization of APC polynomials and matrices
F = zeros(1,na+1);
E = [1 zeros(1,N2-1)];
G = zeros(1,N2+nb-1);
GAMMA   = zeros(N2-nk+1,Nu);
PHI   = zeros(N2,1);
F0    = zeros(N2,1);
rhoI  = rho*eye(Nu);
invHG = zeros(Nu,N2-N1+1);
A     = [1 zeros(1,na)];
B     = zeros(1,nb);

% Miscellaneous initializations
t    = -Ts;
yhat = y_0;

fighandle=progress;

% Initialization of Simulink system
if strcmp(simul,'simulink')
  simoptions = simset('Solver',integrator,'MaxRows',0); % Set integrator opt.
  eval(['[sizes,x0] = ' sim_model '([],[],[],0);']);    % Get initial states
end


%>>>>>>>>>>>>>>>>>    CALCULATE REFERENCE SIGNAL & FILTER IT     <<<<<<<<<<<<<<<<<<
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


%>>>>>>>>>>>   INITIALIZATION OF VECTORS USED FOR STORING PAST DATA    <<<<<<<<<<<<
ref_data    = ref(1:samples);
u_data      = zeros(samples,1);
y_data      = zeros(samples,1);
yhat_data   = zeros(samples,1);
t_data      = zeros(samples,1);
A_data      = zeros(samples,na+1);
B_data      = zeros(samples,nb);


%------------------------------------------------------------------------------
%-------------------         >>>   MAIN LOOP   <<<           ------------------
%------------------------------------------------------------------------------
for i=1:samples,
  t = t + Ts;
  

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
 
  % Predictive controller
  if strcmp(regty,'apc'),
    PHI = PHI+F0*y;         % Complete calculation of the vector PHI
    U=invHG*(ref(i+N1-1:i+N2-nk)-PHI(N1:N2));
    u = U(1)+u_old(1);
        
  % PID controller
  elseif strcmp(regty,'pid'),
    ui = B1*e + I1;
    um = ui;
    if ui<uimin, um=uimin; end
    if ui>uimax, um=uimax; end
    u = (um-I2)*B2 + I2;
    I1 = I1 + (K*e - (ui - um))*A1;
    I2 = I2 + (um - I2)*A2;
  
  % No control
  else
     u = ref(i);
  end
  
  % Make sure control inputs is within limits
  if u>ulim_max,
     u=ulim_max;
  elseif u<ulim_min
     u=ulim_min;
  end
 
  ey = y - yhat;                          % prediction error (a priori)


%>>>>>>>>>>>>>>>>>>>>>>>>>>       TIME UPDATES        <<<<<<<<<<<<<<<<<<<<<<<<<
  y_old    = shift(y_old,y);
  u_old    = shift(u_old,u);
  ref_old  = shift(ref_old,ref(i));


%>>>>>>>>>>>>>>>>>>>       STORE DATA IN DATA VECTORS      <<<<<<<<<<<<<<<<<<<
  u_data(i)       = u;
  y_data(i)       = y;
  yhat_data(i)    = yhat;
  t_data(i)       = t;
  A_data(i,:)     = A;
  B_data(i,:)     = B;  


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


%>>>>>>>>>>>>>>>>>>>   EXTRACT LINEAR PARAMETERS FROM NN   <<<<<<<<<<<<<<<<<<<<
   % Matrix consisting of the partial derivatives of each output with
   % respect to each of the outputs from the hidden neurons
   d21 = W2;
   for j = H_output',
     d21(j,:) = (1-yhat(j)*yhat(j))*W2(j,:);
   end

   % Matrix with partial derivatives of the output from each hidden neurons
   % with respect to each input:
   d10 = W1;
   for j = H_hidden',
     d10(j,:) = (1-y1(j)*y1(j))*W1(j,:);
   end

   % Matrix with partial derivative of each output with respect to each input
   d20 = d21(1:hidden)*d10;

   A(2:na+1) = -d20(1,1:na);
   B = d20(1,na+1:nab);


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
  invHG=inv(H)*GAMMA(N1:N2,:)';


%>>>>>>>>>>>>>>>>>>      WRITE % OF SIMULATION COMPLETED      <<<<<<<<<<<<<<<<<
  progress(fighandle,floor(100*i/samples));
end
%------------------------------------------------------------------------------
%------------------        >>>   END OF MAIN LOOP   <<<       -----------------
%------------------------------------------------------------------------------
