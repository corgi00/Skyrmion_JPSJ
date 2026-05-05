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

    theta = (1 - r / (beta*lambda)) * (pi);
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

ZZ = reshape(LCM, L*L, 1);

%%

xvals = unique(POSITION(:,1));
yvals = unique(POSITION(:,2));

nx = length(xvals);  % 8
ny = length(yvals);  % 8


Zmat = nan(nx, ny);

for i = 1:length(POSITION)
    xi = POSITION(i,1);
    yi = POSITION(i,2);
    zi = ZZ(i);

    ix = find(xvals == xi);
    iy = find(yvals == yi);

    Zmat(ix,iy) = zi;
end


[Xgrid,Ygrid] = meshgrid(xvals, yvals);
scatter(Xgrid(:), Ygrid(:), 60, real(Zmat(:)), 'filled');
axis equal;
axis off;
%% -1

caxis([-1.25, 9]);
cmin = -1.25;
cmax = 9;


n = 256;


frac_neg = abs(cmin) / (abs(cmin) + cmax);  
n_blue = round(n * frac_neg);      
n_red  = n - n_blue;               


r_blue = linspace(0, 1, n_blue)';   
gamma_blue = 0.6;     

t = linspace(0,1,n_blue)'.^gamma_blue;

r_blue = t;
g_blue = t;
b_blue = ones(n_blue,1);

blue = [r_blue, g_blue, b_blue];


target_red = [180, 0, 0] / 255;
r_red = linspace(1, target_red(1), n_red)'; 
g_red = linspace(1, target_red(2), n_red)'; 
b_red = linspace(1, target_red(3), n_red)'; 
red = [r_red, g_red, b_red];


custom_cmap = [blue; red];
colormap(custom_cmap);

cb = colorbar;
cb.Ticks = [-1, 0, 1, 3, 5, 7, 9];
cb.FontSize = 14;

axis off
axis equal


%% ==========

block = 28;                
N = size(Zmat, 1);          


startIdx = floor((N - block)/2) + 1;
endIdx   = startIdx + block - 1;

subZ = Zmat(startIdx:endIdx, startIdx:endIdx);  
subX = xvals(startIdx:endIdx);
subY = yvals(startIdx:endIdx);


[Xsub, Ysub] = meshgrid(subX, subY);

figure;
scatter(Xsub(:), Ysub(:), 60, real(subZ(:)), 'filled');
%%

colormap(blue);

caxis([cmin 0]);   
colorbar;
cb2 = colorbar;
cb2.FontSize = 14;     

axis equal;
axis off;


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

