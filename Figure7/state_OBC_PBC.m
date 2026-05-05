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

All_BOTT = [];   % 假设 CELL_N 从 1 到 16


CELL_N=9;
CELL=kron(eye(CELL_N),h1)+kron(diag(ones(CELL_N-1,1),1),t1)+kron(diag(ones(CELL_N-1,1),-1),t1');
CELL_t2=kron(eye(CELL_N),t2);
CELL_Z=zeros(size(CELL));

H=kron(eye(CELL_N),CELL)+kron(diag(ones(CELL_N-1,1),1),CELL_t2)+...
    kron(diag(ones(CELL_N-1,1),-1),CELL_t2');

%%
[Psi,E1]=eig(H);
BOTT_state=[];

E=diag(E1);
E_obc=E;

halfN=150;
% ===== OBC（红色）=====
h_obc_scatter = scatter(1:halfN, E_obc(1:halfN), 0.5, 'r', 'filled');
hold on
h_obc_line = plot(1:halfN, E_obc(1:halfN), 'r');

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
CELL_N=9;
CELL=kron(eye(CELL_N),h1)+kron(diag(ones(CELL_N-1,1),1),t1)+kron(diag(ones(CELL_N-1,1),-1),t1');
CELL_t2=kron(eye(CELL_N),t2);
CELL_Z=zeros(size(CELL));

H=kron(eye(CELL_N),CELL)+kron(diag(ones(CELL_N-1,1),1),CELL_t2)+...
    kron(diag(ones(CELL_N-1,1),-1),CELL_t2');
%% gama

dt_=[1 0;-1  0;0 1; 0 -1];
kt_=zeros(4,2);
intra=H;
LE=(length(H));

AC3=kron(eye(4),hop);
tm1_CELL=zeros(4);
tm1_CELL(4,1) = 1;
tm1_cell=kron(tm1_CELL,AC3);
A = zeros(CELL_N,CELL_N);
A(CELL_N,1) = 1;
tm1=kron(kron(A,diag(ones(CELL_N,1))),tm1_cell);
tm2=tm1';

tm3_cell=kron(eye(4),kron(tm1_CELL,hop));
tm3=kron(kron(eye(CELL_N),A),tm3_cell);
tm4=tm3';

TM1={tm1};TM2={tm2};TM3={tm3};
TM4={tm4};
TM=[TM1;TM2;TM3;TM4];

bloch=[];
for i=1:1:length(kt_(:,1))
dt=(dt_(i,:));kt=(kt_(i,:));
result = exp(1i * dot(kt, dt) * pi * 2);
bloch=[bloch,result];
end

tk_=[];
for i=1:1:length(kt_(:,1))
TK=bloch(:,i)*TM{i,:};
tk_=[tk_;TK];
end

rows_per_cell = LE;
num_cells = size(tk_, 1) / rows_per_cell;  
tk = mat2cell(tk_, repmat(rows_per_cell, num_cells, 1), size(tk_, 2));
group_size = 4;
num_groups = num_cells / group_size;
final_tk = cell(num_groups, 1);
for i = 1:num_groups
    start_idx = (i-1)*group_size + 1;
    end_idx = i*group_size;

    sum_matrix = zeros(size(tk{start_idx}));

    for j = start_idx:end_idx
        sum_matrix = sum_matrix + tk{j};
    end

    final_tk{i} = sum_matrix;
end
resultCells = final_tk;

RC_=[];
for i=1:1:length(resultCells)
RC=intra+resultCells{i};
RC_=[RC_;RC];
end
%%

[Psi,E1]=eig(RC_);
E=diag(E1);
E_pbc=E;
N = length(E);
% ===== PBC（蓝色）=====
h_pbc_scatter = scatter(1:halfN, E_pbc(1:halfN), 0.5, 'b', 'filled');
h_pbc_line = plot(1:halfN, E_pbc(1:halfN), 'b');
%%
xlim([0 halfN])
xticks([0 30 60 90 120 halfN])
%%
y1 = -12.8856;
y2 = -12.6683;

xlim_vals = xlim;   % 取得当前 x 轴范围

patch([xlim_vals(1) xlim_vals(2) xlim_vals(2) xlim_vals(1)], ...
      [y1 y1 y2 y2], ...
      [0.8 0.8 0.8], ...      % 灰色
      'EdgeColor','none', ...
      'FaceAlpha',0.5);       % 透明度

hold on
%%
set(gca, 'FontWeight', 'bold', 'FontSize', 14);
% ========= 只框住坐标轴绘图区域（不含轴外文字）
box on
set(gca, 'LineWidth', 1.2)
%%
% ===== legend 只绑定“线”更好看 =====
lgd = legend([h_pbc_line, h_obc_line], {'PBC', 'OBC'}, ...
    'Location', 'southeast');

set(h_pbc_line, 'Marker','o', 'MarkerSize',4, ...
    'MarkerFaceColor','b')

set(h_obc_line, 'Marker','o', 'MarkerSize',4, ...
    'MarkerFaceColor','r')

lgd.Box = 'on';
lgd.Color = [1 1 1 0.85];
lgd.EdgeColor = [0.8 0.8 0.8];
lgd.FontSize = 18;
lgd.FontWeight = 'normal';   % ← 关键这一行
%%

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

