clear all
close all
clc

%% Chose demo
flagGraph = 2; % 1 Unsigned graph；2 Signed graph
flagTest = 2;  % 1 Synchronous；  2 Asynchronous

%% Parameters
delta = 0.45;

%% initial value
iter = 51;
n = 8;
p = 3;
G_A = zeros(n,n);
G_A(1,5) = -1; G_A(1,8) = 1;
G_A(2,1) = 1;  G_A(2,8) = -1;
G_A(3,2) = -1; G_A(3,8) = 1;
G_A(4,3) = 1; G_A(4,5) = -1;
G_A(5,3) = 1; G_A(5,6) = -1;
G_A(6,1) = 1; G_A(6,5) = -1;
G_A(7,1) = 1; 
L = diag(sum(G_A,2)) - G_A;
lamda = eig(L);

switch flagGraph
    case 1
    G_A = abs(G_A); % no-signed
    case 2
    G_A = (G_A); % signed
end


%% Define neighbors
GA1 = G_A(1:n,1:n);  
GA1((GA1<0)) = 0;
DA1 = diag(sum(GA1,2)); 
GA2 = G_A(1:n,1:n);  
GA2((GA2>0)) = 0;
DA2 = diag(sum(GA2,2));  

%% Asynchronous setting
array_length = iter;
max_interval = 2;  % Asynchronous
array = zeros(n, array_length); 
array(1:n,1) = 1; 
for i = 1 : n
current_index = 1;
while current_index < array_length
    next_index = current_index + randi(max_interval);
    if next_index <= array_length
        array(i,next_index) = 1;
        current_index = next_index;
    else
        break; 
    end
end
end

%% System Model
A = [1.1 -0.2 0.2; 0.5 0.7 -0.3; -0.4 0.4 1];
% A = [0.75 -0.2 0.2; 0.5 0.7 -0.3; -0.4 0.4 1]; % case 2
B = [1 0 0 1; 0 1 0 0; 0 0 1 0];
K = delta*B'*((B*B')^(-1))*A;

%% Leader and Followers
XL = zeros(p,iter);
XF = zeros(p*(n-1),iter);
UF = zeros(p*(n-1),iter); 
E = zeros(p*(n-1),iter);
XL(:,1) = 10*rand(p,1)-5;
XF(:,1) = 10*rand(1,p*(n-1))-5;
E(:,1) = XF(:,1) - kron(ones(n-1,1),XL(:,1));
X = zeros(p*(n),iter);
X(:,1) = [XF(:,1);XL(:,1)];
U = zeros(p*(n),iter); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
G_A_original = G_A;

for k = 1:(iter-1)
    
F = zeros(n,n);
for i = 1 : n
    for j = 1 : n
   if G_A(i,j)>0
           Ef = abs(X(1+p*(j-1),k) - X(1+p*(i-1),k));
           F(i,j) = 0.05 * 2^(-Ef) + 0.95;
   else if G_A(i,j)<0
           Ef = abs(X(1+p*(j-1),k) - X(1+p*(i-1),k));
           F(i,j) = 0.02/(1+Ef) + 0.06;    
       else
           F(i,j) = 0;
       end
   end
    end
end

G_A = (G_A_original) .* F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = diag(sum(G_A,2)) - G_A;

switch flagTest
    case 1 % Synchronous
%% 
H1 = eye(n-1) - delta*L(1:n-1,1:n-1);
H = [H1,-delta*L(1:n-1,n); zeros(1,n-1) 1];
X(:,k+1) = kron(H, A) * X(:,k);
U(:,k) = -kron(L, B*K) * X(:,k);

    case 2 % Asynchronous
%% 
H1 = eye(n-1) - delta*L(1:n-1,1:n-1).*array(1:n-1,k);
H = [H1,-delta*L(1:n-1,n).*array(1:n-1,k); zeros(1,n-1) 1];
X(:,k+1) = kron(H, A) * X(:,k);
U(:,k) = -kron(L, B*K) .* kron(ones(3,1),array(1:n,k)) * X(:,k);
end

%% Error
XF(:,k) = X(1:p*(n-1),k);
XL(:,k) = X(p*(n-1)+1:p*(n),k);
E(:,k) = XF(:,k) - kron(ones(n-1,1),XL(:,k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
%% Position
for i = 1 : n
plot(0:1:iter-1, X(1+3*(i-1),:),'linewidth',1.5); hold on
end
ylabel('The trajectories of agents $x_{i1}$', 'Interpreter', 'latex','FontSize',14);
xlabel('Time/$t$', 'Interpreter', 'latex','FontSize',14);
legend('Follower1','Follower2','Follower3','Follower4','Follower5','Follower6','Follower7','Leader', 'Interpreter', 'latex','FontSize',10);
set(gca,'FontName','Times New Roman');
grid on

figure(2)
%% Velocity
for i = 1 : n
plot(0:1:iter-1, X(2+3*(i-1),:),'linewidth',1.5); hold on
end
ylabel('The trajectories of agents $x_{i2}$', 'Interpreter', 'latex','FontSize',14);
xlabel('Time/$t$', 'Interpreter', 'latex','FontSize',14);
legend('Follower1','Follower2','Follower3','Follower4','Follower5','Follower6','Follower7','Leader', 'Interpreter', 'latex','FontSize',10);
set(gca,'FontName','Times New Roman');
grid on

figure(3)
%% Accelrate
for i = 1 : n
plot(0:1:(iter-1), X(3+3*(i-1),:),'linewidth',1.5); hold on
end
ylabel('The trajectories of agents $x_{i3}$', 'Interpreter', 'latex','FontSize',14);
xlabel('Time/$t$', 'Interpreter', 'latex','FontSize',14);
legend('Follower1','Follower2','Follower3','Follower4','Follower5','Follower6','Follower7','Leader', 'Interpreter', 'latex','FontSize',10);
set(gca,'FontName','Times New Roman');
grid on

figure(4)
%% Error
for i = 1 : n-1
plot(0:1:(iter-1), E(1+3*(i-1),:),'linewidth',1.5); hold on
end
ylabel('The trajectories of agents $x_{i3}$', 'Interpreter', 'latex','FontSize',14);
xlabel('Time/$t$', 'Interpreter', 'latex','FontSize',14);
legend('Follower1','Follower2','Follower3','Follower4','Follower5','Follower6','Follower7', 'Interpreter', 'latex','FontSize',10);
set(gca,'FontName','Times New Roman');
grid on

figure(5)
%% Error
for i = 1 : n-1
plot(0:1:(iter-1), E(2+3*(i-1),:),'linewidth',1.5); hold on
end
ylabel('The trajectories of agents $x_{i3}$', 'Interpreter', 'latex','FontSize',14);
xlabel('Time/$t$', 'Interpreter', 'latex','FontSize',14);
legend('Follower1','Follower2','Follower3','Follower4','Follower5','Follower6','Follower7', 'Interpreter', 'latex','FontSize',10);
set(gca,'FontName','Times New Roman');
grid on

figure(6)
%% Error
for i = 1 : n-1
plot(0:1:(iter-1), E(3+3*(i-1),:),'linewidth',1.5); hold on
end
ylabel('The trajectories of agents $x_{i3}$', 'Interpreter', 'latex','FontSize',14);
xlabel('Time/$t$', 'Interpreter', 'latex','FontSize',14);
legend('Follower1','Follower2','Follower3','Follower4','Follower5','Follower6','Follower7', 'Interpreter', 'latex','FontSize',10);
set(gca,'FontName','Times New Roman');
grid on