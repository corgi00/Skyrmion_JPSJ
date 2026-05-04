clc
clear

%% 
DT=readmatrix('dt.txt');
DT_=(reshape(DT.',2,16024/2))';
dt_=[DT_];

basic=[1 0; -1 0; 0 1; 0 -1];
repeated_basic = repmat(basic, 400, 1);
dt_=repeated_basic;
%%
KT=readmatrix('kt_whole.txt');
kt_=(reshape(KT',2,length(KT(:,1))/2))';
%%

tm=readmatrix('tm_n4.txt');
 rows_per_cell = 32;
num_cells = size(tm, 1) / rows_per_cell; 

TM_ = mat2cell(tm, repmat(rows_per_cell, num_cells, 1), size(tm, 2));
TM=[TM_];
%%

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

%% BERRY_CURVATURE

Chern=[];
for n=1:1:16

f=[];

PHI=[];
for i=4:4:400

% n=1:16;

[v1,e1] = eig(H{i-3,1});
wf1 = v1(:, n)';

[v2,e2] = eig(H{i-2,1});
wf2 = v2(:, n)';

[v3,e3] = eig(H{i-1,1});
wf3 = v3(:, n)';

[v4,e4] = eig(H{i,1});
wf4 = v4(:, n)';


out1 = conj(wf1) * wf2.';
out2 = conj(wf2) * wf3.';
out3 = conj(wf3) * wf4.';
out4 = conj(wf4) * wf1.';

m=out1 * out2 * out3 * out4;
d = det(m); 

dr=real(d);
di=imag(d);
phi = atan2(di,dr)./(4*0.05*0.05);
PHI=[PHI,phi];

end

S=(sum(PHI))./(2*pi*100);

Chern=[Chern,S];

end
C_sum=sum(Chern)


scatter(Chern, 1:numel(Chern), 'filled');
hold on
plot(Chern, 1:numel(Chern));
