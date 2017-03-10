%Î°´¨
Aforlin = [1 2];
A = [1 Aforlin];
B = [3 4];
rhoI = 0.1;
a1 = Aforlin(1);
a2 = Aforlin(2);
b0 = B(1);
b1 = B(2);
u_old = [1 2]';
y = 1;
y_old = [3 4]';
difu = 0;
ref = ones(7,1);
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
M = G'*G+rhoI;
invMG = inv(M)*G'
U=invMG*(ref-S*[y;y_old(1);y_old(2)]-H*difu)