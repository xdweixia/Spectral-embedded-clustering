% run SEC_YTF
% ==================load_data===============================
clc;
clear all;
X = load('YTF_X.mat');     % X: each column is a data point
Label=load('YTF_Label');
X = X.X;
Y = Label.L + 1;

[m,n]=size(X);
result=zeros(n,n);
xi = 5;
for i=1:n
	for j=1:n
		a=X(:,i)-X(:,j);
		result(i,j)=norm(a,'fro');
		result(i,j) = exp(-result(i,j)/xi);
		if i == j
            result(i,j) = 0;
        end
	end
end 
% D=diag(sum(result));
% L=D-result;
N = size(result,1);
DN = diag( 1./sqrt(sum(result)+eps) );
L = speye(N) - DN * result * DN;
N_cluster = 41;
gamma = 1;
% -10,-7,-4,-1,2,5,8
mu = 1e+5;
para.gamma=gamma;
para.mu=mu;
para.nc=N_cluster;
[F, W, b] = SEC(X, L, para);
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans
groups = kmeans(F,N_cluster,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
[A, nmi, avgent] = compute_nmi (Y, groups);
ACC = Accuracy(Y, groups);