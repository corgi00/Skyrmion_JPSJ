clc
clear

DT=readmatrix('dt.txt');
dt_=(reshape(DT.',2,16024/2))';


gama=[0,0]; X=[0.5,0]; M=[0.5,0.5];

X1 = linspace(gama(1,1), X(1,1), 100);
Y1 = linspace(gama(1,2), X(1,2), 100);

X2 = linspace(X(1,1), M(1,1), 100);
Y2 = linspace(X(1,2), M(1,2), 100);

X3 = linspace(M(1,1), gama(1,1), 100);
Y3 = linspace(M(1,2), gama(1,2), 100);

KT=[X1 X2 X3; Y1 Y2 Y3];

kt_=repelem(KT', 4, 1);

tm=readmatrix('tm_N4.txt');
rows_per_cell = 32;
num_cells = size(tm, 1) / rows_per_cell;  
TM = mat2cell(tm, repmat(rows_per_cell, num_cells, 1), size(tm, 2));

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

num_cells = size(tk_, 1) / rows_per_cell;  
tk=mat2cell(tk_, repmat(rows_per_cell, num_cells, 1), size(tk_, 2));


originalCells = tk;


resultCells = cell(ceil(numel(originalCells)/4), 1);


for i = 1:4:numel(originalCells)
   
    sumMatrix = zeros(32, 32);

    
    endIndex = min(i+3, numel(originalCells));

    
    for j = i:endIndex
        
        sumMatrix = sumMatrix + originalCells{j};
    end

  
    resultCells{ceil((i-1)/4) + 1} = sumMatrix;
end
%%

X=[];
Y=[];
LATTICE=1.5;
for i=-LATTICE:1:LATTICE
    for j= -LATTICE:1:LATTICE
X=[X,i];
Y=[Y,j];
    end
end
position=[X;Y]';


t=1;m=1;gamma=0;beta=1;
hop=[t 0;0 t];


J=-10;

sx = [0 1; 1 0];
            
sy = [0 -1.0i; 1.0i 0];
      
sz = [1 0; 0 -1];



M=[];
N=length(position(:,1));

for x=1:1:N

    lamda=round(LATTICE);
    X=position(x,1);Y=position(x,2);

    r  = sqrt(X^2 + Y^2);
    theta = (1 - r / (beta*lamda)) * pi;
    theta = max(min(theta, pi), 0);

    phi_xy = m*atan2(Y,X)+gamma;

    nr(1) = sin(theta) * cos(phi_xy);
    nr(2) = sin(theta) * sin(phi_xy);
    nr(3) = cos(theta); 
    nr = nr * J;
    n=sx*nr(:,1) + sy*nr(:,2) + sz*nr(:,3);

M=[M,n];%

    end


A = M;  


D = zeros(N*2, N*2);


for k = 1:N
 
    cols = A(:, 2*k-1:2*k);


    D(2*k-1:2*k, 2*k-1:2*k) = cols;
end


HOP1=kron(diag((ones(sqrt(N)-1, 1)),-1),hop)+kron(diag((ones(sqrt(N)-1, 1)),1),hop);
HOP1_=kron(diag((ones(sqrt(N), 1))),HOP1);

HOP2=kron(diag((ones(sqrt(N), 1))),hop);
HOP2_=kron(diag((ones(sqrt(N)-1, 1)),-1),HOP2)+kron(diag((ones(sqrt(N)-1, 1)),1),HOP2);

HH=D+HOP1_+HOP2_;


intra=HH;
%%
RC_=[];
for i=1:1:length(resultCells)
RC=intra+resultCells{i};
RC_=[RC_;RC];
end

mout=RC_;
num_cells = size(mout, 1) / rows_per_cell;  
H = mat2cell(mout, repmat(rows_per_cell, num_cells, 1), size(mout, 2));

E=[];
for i=1:1:300
[w,e]=eig(H{i,:});
e_=diag(e);
E=[E;e_];
end


num_cells = size(E, 1) / rows_per_cell;  
EE = mat2cell(E, repmat(rows_per_cell, num_cells, 1), size(E, 2));

EE_ = cell(300, 1); 
for i = 1:1:300
    test = EE{i, 1}; 
    test_ = test(1:16, :); 
    EE_{i, 1} = test_;
end


x=[];
for i=1:1:300
x=[x,i];
x_=x';
end
X=(repmat(x_, 1, 16))';

for i=1:1:300
scatter(X(:,i),EE_{i,:},2,'blue')
hold on
end 

xticks([0 100 200 300]);
xticklabels({'\Gamma' 'X' 'M' '\Gamma'});
set(gca, 'FontWeight', 'bold', 'FontSize', 14);

set(gcf,'Color','w');


%%

ax = gca;
pos = ax.Position;   

annotation('rectangle', ...
    pos, ...
    'LineWidth', 1.2, ...
    'Color', 'k');
