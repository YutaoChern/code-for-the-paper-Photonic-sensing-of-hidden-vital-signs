function [p,s]=MYIAA(x,K) %MIAA：多快拍IAA（是指对多列数据进行IAA估计后取平均，来获得更稳定的估计结果）
[M,N]=size(x);%注意 要求输入格式为M=数据长度，N为快拍数
if nargin==1
    K=M*1;
    CC=10;K=K*CC;
end
t=(0:M-1)';

w=(0:K-1)'/K*2*pi;
A=exp(1j*t*w.');
% A0=A(:,16*CC:26*CC);A0=[A(:,1:CC:16*CC) A0 A(:,26*CC:CC:end)];A0=A;K=length(A0);
x0=fft(x,K)/M;
p=sum(abs(x0).^2,2)/N;p=ones(K,1);
Maxiter=100;
z=zeros(size(x));
s=zeros(size(x));
for i=1:Maxiter
R = (A .* p')*A'+0*x*x';
Ri=inv(R);
 X=A'*Ri;  
% RAA = X * A; 
diag_C = sum(X .* A.', 2); 
S = (X*x) ./ diag_C; 
p = sum(abs(S).^2 ,2)/ N; % K x 1 向量
p=p/max(p);
% toc
    if abs(norm(s)-norm(z))/norm(s)<0.001
        break;
    else
        z=s;
    end
end


