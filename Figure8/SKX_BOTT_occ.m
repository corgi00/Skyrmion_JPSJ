clc
clear

%% 1
%% ——— ① 生成 0~3 的 4×4 网格 ———
xlist = 0.5 : 3.5;              % [-4, -3, -2, -1] 共 4 点
ylist = 0.5 : 3.5;  
[Xg, Yg] = meshgrid(xlist, ylist);
position1 = [Xg(:), Yg(:), zeros(numel(Xg),1)];   % 16 × 3

%% ——— ② 设 Skyrmion 半径和中心 ———
lambda = 2;                 % 半径调小些，刚好落在 4×4 网格里
beta   = 1;
x0 = 2;                   % 中心放在网格几何中心 (1.5,1.5)
y0 = 2;
m  = 1;
gamma = 0;
J=-10; 
t=1;
hop=[t 0;0 t];
sx = [0 1; 1 0];
            
sy = [0 -1.0i; 1.0i 0];
      
sz = [1 0; 0 -1];

%% ——— ③ 计算自旋方向 ———
M=[];

for idx = 1:size(position1,1)
    Xp = position1(idx,1);
    Yp = position1(idx,2);

    dx = Xp - x0;
    dy = Yp - y0;
    r  = sqrt(dx^2 + dy^2);

    theta = (1 - r / (beta*lambda)) * pi;
    theta = max(min(theta, pi), 0);      % 截断到 [0,π]
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
Z=zeros(size(h1));
M=[1,5,9,13];
All_BOTT = [];   % 假设 CELL_N 从 1 到 16
for m=1:1:4

CELL_N=M(m);
CELL=kron(eye(CELL_N),h1)+kron(diag(ones(CELL_N-1,1),1),t1)+kron(diag(ones(CELL_N-1,1),-1),t1');
CELL_t2=kron(eye(CELL_N),t2);
CELL_Z=zeros(size(CELL));

H=kron(eye(CELL_N),CELL)+kron(diag(ones(CELL_N-1,1),1),CELL_t2)+...
    kron(diag(ones(CELL_N-1,1),-1),CELL_t2');

%%
[Psi,E1]=eig(H);
BOTT_OCC=[];

for occ=1:1:16
E=diag(E1);
%Psi_=Psi(:,1:100);
Psi_=Psi(:,1:CELL_N*CELL_N*occ);
P=Psi_*Psi_';% P
Q=eye(size(P))-P; %Q
%%
cellL = 4;    % 每个原胞 4×4
nx    = CELL_N;    % x 方向 4 个原胞
ny    = CELL_N;    % y 方向 4 个原胞

POSITION = generate_supercell_ordered(cellL, nx, ny);

POSITION_=repelem(POSITION,2,1);
X_res=POSITION_(:,1);
Y_res=POSITION_(:,2);

eps_norm = 1e-6;                           
X_res = (1-eps_norm) * ...
        (X_res - min(X_res)) / (max(X_res) - min(X_res));
Y_res = (1-eps_norm) * ...
        (Y_res - min(Y_res)) / (max(Y_res) - min(Y_res));
%%
U=P*(diag(exp(1i*2*pi*X_res)))*P+Q;
V=P*(diag(exp(1i*2*pi*Y_res)))*P+Q;

se=eig(V*U*V'*U');
Bott=imag(sum(log(se)))/(2*pi);
BOTT_OCC=[BOTT_OCC,Bott];
end
All_BOTT = [All_BOTT;BOTT_OCC];

end
figure('Position',[100 100 800 400]);

colors = lines(4);  
hold on

order = [3,2,4,1];
hBott = gobjects(1,4);  % Bott index 句柄

for i = 1:4
    if i == 3
        thisColor = [0 0.6 0];   % 绿色
    else
        thisColor = colors(order(i),:);
    end
    
    hBott(i) = plot(1:length(All_BOTT(i,:)), All_BOTT(i,:), ...
                    'Color', thisColor, 'LineWidth', 2,'LineStyle', '--');
    scatter(1:length(All_BOTT(i,:)), All_BOTT(i,:), ...
            'filled', 'MarkerFaceColor', thisColor, ...
            'MarkerEdgeColor', thisColor);
end

function POSITION = generate_supercell_ordered(cellL, nx, ny)

block_pts = cellL^2;                 % 每个原胞的点数 (16)

POSITION = zeros(nx*ny*block_pts, 3); % 预分配
ptr = 1;                              % 写入指针

for ix = 0:nx-1          % 外层先扫 x-方向的大原胞
    xlist = ix*cellL + (0.5 : 1 : cellL-0.5);
    for iy = 0:ny-1      % 内层再扫 y-方向的大原胞
        ylist = iy*cellL + (0.5 : 1 : cellL-0.5);
        [Xg, Yg] = meshgrid(xlist, ylist);
        blk      = [Xg(:), Yg(:), zeros(block_pts,1)];
        POSITION(ptr : ptr+block_pts-1, :) = blk;
        ptr = ptr + block_pts;
    end
end
end

