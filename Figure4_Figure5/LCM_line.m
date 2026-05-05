clc
clear

%% ——————
xlist = 0.5 : 3.5;             
ylist = 0.5 : 3.5;  
[Xg, Yg] = meshgrid(xlist, ylist);
position1 = [Xg(:), Yg(:), zeros(numel(Xg),1)];   

%% ——————
lambda = 2;             
beta   = 1;
x0 = 2;              
y0 = 2;
m  = 1;
gamma = 0;
J=-10; 
t=1;
hop=[t 0;0 t];
sx = [0 1; 1 0];
            
sy = [0 -1.0i; 1.0i 0];
      
sz = [1 0; 0 -1];

%% ——————
M=[];

for idx = 1:size(position1,1)
    Xp = position1(idx,1);
    Yp = position1(idx,2);

    dx = Xp - x0;
    dy = Yp - y0;
    r  = sqrt(dx^2 + dy^2);

    theta = (1 - r / (beta*lambda)) * pi;
    theta = max(min(theta, pi), 0);      
    phi   = m * atan2(dy, dx) + gamma;

    nr = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    nr = nr / norm(nr);

    nr(1) = sin(theta) * cos(phi);
    nr(2) = sin(theta) * sin(phi);
    nr(3) = cos(theta); 
    nr = nr * J;
    n=sx*nr(:,1) + sy*nr(:,2) + sz*nr(:,3);

M=[M,n];%
end

A = M; 
N=length(position1(:,1));
D = zeros(N*2, N*2);

for k = 1:N
    cols = A(:, 2*k-1:2*k);

    D(2*k-1:2*k, 2*k-1:2*k) = cols;
end

HOP1=kron(diag((ones(sqrt(N)-1, 1)),-1),hop)+kron(diag((ones(sqrt(N)-1, 1)),1),hop);
HOP1_=kron(diag((ones(sqrt(N), 1))),HOP1);

HOP2=kron(diag((ones(sqrt(N), 1))),hop);
HOP2_=kron(diag((ones(sqrt(N)-1, 1)),-1),HOP2)+kron(diag((ones(sqrt(N)-1, 1)),1),HOP2);

h1=D+HOP1_+HOP2_;
%%
T1=zeros(4,4);
T1(4,1) = 1;
t1=kron(eye(4),(kron(T1,hop)));

T2=zeros(4);
T2(4,1) = 1;
t2=kron(T2,(kron(eye(4),hop)));
%%
AV=[];

Z=zeros(size(h1));
CELL_N=11;
CELL=kron(eye(CELL_N),h1)+kron(diag(ones(CELL_N-1,1),1),t1)+kron(diag(ones(CELL_N-1,1),-1),t1');
CELL_t2=kron(eye(CELL_N),t2);
CELL_Z=zeros(size(CELL));

H=kron(eye(CELL_N),CELL)+kron(diag(ones(CELL_N-1,1),1),CELL_t2)+...
    kron(diag(ones(CELL_N-1,1),-1),CELL_t2');
%%
L=sqrt(length(H)./2);
[Psi,E1]=eig(H);

E=diag(E1);

Psi_=Psi(:,1:CELL_N*CELL_N);
P=Psi_*Psi_';% P
Q=eye(size(P))-P; %Q
%%
cellL = 4;    
nx    = CELL_N;    
ny    = CELL_N;   

POSITION = generate_supercell_ordered(cellL, nx, ny);

POSITION_=repelem(POSITION,2,1);
X0=(POSITION_(:,1));
Y0=(POSITION_(:,2));
X=diag(X0);
Y=diag(Y0);

PxP=P*X*P; PyP=P*Y*P;
commutator=PxP * PyP - PyP * PxP;
lcm=diag(-2 * pi * 1i*commutator);
LCM=sum(reshape(lcm,2,L*L));

%% reshape_LCM
A = LCM;
n_split = CELL_N*CELL_N;
chunk_size = 16;

B = cell(1, n_split); 

for i = 1:n_split
    idx_start = (i-1)*chunk_size + 1;
    idx_end = i*chunk_size;
    B{i} = A(idx_start:idx_end);
end

block_size = 4;  


M = zeros(CELL_N * block_size, CELL_N * block_size);

for i = 1:CELL_N
    for j = 1:CELL_N
        idx = (i - 1) * CELL_N + j;
        subblock = reshape(B{idx}, block_size, block_size)';  
        row_range = (i-1)*block_size + (1:block_size);
        col_range = (j-1)*block_size + (1:block_size);
        M(row_range, col_range) = subblock;
    end
end

%% reshape_POSITION_X
POSITION_X=POSITION(:,1);
A = POSITION_X;
n_split = CELL_N*CELL_N;
chunk_size = 16;

MX = cell(1, n_split);  

for i = 1:n_split
    idx_start = (i-1)*chunk_size + 1;
    idx_end = i*chunk_size;
    MX{i} = A(idx_start:idx_end);
end

block_size = 4;  
MMX = zeros(CELL_N * block_size, CELL_N * block_size);

for i = 1:CELL_N
    for j = 1:CELL_N
        idx = (i - 1) * CELL_N + j;
        subblock = reshape(MX{idx}, block_size, block_size)';  
        row_range = (i-1)*block_size + (1:block_size);
        col_range = (j-1)*block_size + (1:block_size);
        MMX(row_range, col_range) = subblock;
    end
end

%% reshape_POSITION_Y
POSITION_Y=POSITION(:,2);
A = POSITION_Y;
n_split = CELL_N*CELL_N;
chunk_size = 16;

MY = cell(1, n_split);  

for i = 1:n_split
    idx_start = (i-1)*chunk_size + 1;
    idx_end = i*chunk_size;
    MY{i} = A(idx_start:idx_end);
end

MMY = zeros(CELL_N * block_size, CELL_N * block_size);

for i = 1:CELL_N
    for j = 1:CELL_N
        idx = (i - 1) * CELL_N + j;
        subblock = reshape(MY{idx}, block_size, block_size)'; 
        row_range = (i-1)*block_size + (1:block_size);
        col_range = (j-1)*block_size + (1:block_size);
        MMY(row_range, col_range) = subblock;
    end
end
%% plot

row_idx = 23;%line
x_data = MMX(:,1) + 0.5;
y_data = M(row_idx,:);


colors = repmat([0 0 1], length(y_data), 1); 
colors(y_data > 0, :) = repmat([1 0 0], sum(y_data > 0), 1);  %

%% plot
scatter(x_data, y_data, 40, colors, 'filled');
hold on
plot(x_data, y_data, 'k');  
figure(1)
pbaspect([6 2 1])
ylim([-1.6 9]);
xlim([0 45]);%原图
set(gca, 'FontSize', 12);   
yticks([-1 0 1 3 5 7 9]);
xticks([0 5 10 15 20 25 30 35 40 45]);


figure(2)
plot(x_data, y_data, 'k');  
hold on
scatter(x_data, y_data, 40, colors, 'filled');
pbaspect([6 2 1])
ylim([-1.3 -0.6]);
xlim([16.5 28.5]);%放大
set(gca, 'FontSize', 16);   
box on
%%
function POSITION = generate_supercell_ordered(cellL, nx, ny)
block_pts = cellL^2;                 

POSITION = zeros(nx*ny*block_pts, 3); 
ptr = 1;                              

for ix = 0:nx-1         
    xlist = ix*cellL + (0.5 : 1 : cellL-0.5);
    for iy = 0:ny-1     
        ylist = iy*cellL + (0.5 : 1 : cellL-0.5);
        [Xg, Yg] = meshgrid(xlist, ylist);
        blk      = [Xg(:), Yg(:), zeros(block_pts,1)];
        POSITION(ptr : ptr+block_pts-1, :) = blk;
        ptr = ptr + block_pts;
    end
end
end

