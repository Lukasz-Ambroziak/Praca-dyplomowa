function [A] = thinning(K)
% K-j¹dro brzegowe procesu
% funkcja zwraca zbiór A, wylosowany zgodnie z j¹drem brzegowym K
[m,n]=size(K);
if m~=n
    disp("z³y rozmiar macierzy");
end
[T,w]=chol(eye(n)-K);
if w>0
    P=(w:n);
else
    P=[];
    w=n+1;
end
X=zeros(1,n);
x=0;
T=T';
O=inv(T);
q=ones(n,1);
z=zeros(n,1);
for i=1:w-1%próbkowanie procesu dominuj¹cego
    p=rand;
    z=O(1:i-1,1:i-1)*K(1:i-1,i);
    q(i)=K(i,i)+z'*z;
    if p<q(i) 
        x=x+1;
        X(1,x)=i;
    end
end
p=size(P,2);
if p>0
    X(1,x+1:x+p)=P;
    x=x+p;
end
A=zeros(1,n);
a=0;
if x==0
    A=[];
    return;
end
B=zeros(1,n);
if X(1,1)>1
    B(1,1:X(1)-1)=(1:X(1)-1);
end
C=[];
v=X(1,1);
czypier=1;
T=zeros(n);
t=X(1,1)-1;
D=zeros(1,n);
for j=X(1,1:x)    
    if czypier==0
        C=(v+1:j-1);   
        if t>0 && size(C,2)>0 %aktualizacja dekompozycji choleskiego
            o=size(C,2);
            T(t+1:t+o,1:t)=(T(1:t,1:t)\K(B(1,1:t),C))';%%%%
            T(t+1:t+o,t+1:t+o)=(chol(K(C,C)-K(C,B(1,1:t))*((eye(t)-K(B(1,1:t),B(1,1:t)))\K(B(1,1:t),C))))';
            t=t+o;
        end
    else
        T(1:t,1:t)=(chol(eye(t)-K(B(1,1:t),B(1,1:t))))';  
        czypier=0;
    end
    if t==0 && size(C,2)>0
        t=size(C,2);
        T(1:t,1:t)=(chol(eye(size(C,2))-K(C,C)))';
    end
    if size(C,2)>0  
        B(1,t-size(C,2)+1:t)=C;
    end
    D(1,a+1)=j;
    if t>0
        J=T(1:t,1:t)\K(B(1,1:t),D(1,1:a+1));
        H=K(D(1,1:a+1),D(1,1:a+1))+J'*J;
    else
        H=K(D(1,1:a+1),D(1,1:a+1));
    end
    p=H(a+1,a+1)-H(a+1,1:a)*(H(1:a,1:a)\H(1:a,a+1));
    r=rand;
    if r<(p/q(j))% wybieramy, które elementy ze zbioru X trafiaj¹ do A, a które do B
        a=a+1;
        A(1,a)=j;
    else
        if t==0
            t=1;
            T(1,1)=sqrt(1-K(j,j));  
        else
            t=t+1;
            T(t,t)=sqrt(1-K(j,j));
            T(t,1:t-1)=(T(1:t-1,1:t-1)\K(B(1,1:t-1),j))';
        end
        B(1,t)=j;
    end
    v=j;
end
A=A(1,1:a);
end




