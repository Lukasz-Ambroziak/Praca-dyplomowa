function [A] = naive(K)
%funkcja zwraca losowy zbiór z procesu wyznaznikowego przy danym j¹drze brzegowym
[m,n]=size(K);
if m~=n
    disp("z³y rozmiar macierzy");
end
A=zeros(1,n);
a=0;
B=zeros(1,n);
b=0;
T=zeros(n);
t=0;
J=[];
for i=1:n % pêtla która wybiera kolejne punkty 
    D=union(A(1,1:a),i);
    if t>0
        J=T(1:t,1:t)\K(B(1,1:b),D);
        H=K(D,D)+J'*J;
    else
        H=K(D,D);
    end
    p=rand;
    if p<H(a+1,a+1)-H(a+1,1:a)*(H(1:a,1:a)\H(1:a,a+1))
        a=a+1;
        A(1,a)=i;
    else
        if t==0
            t=1;
            T(t,t)=sqrt(1-K(i,i));
        else
            t=t+1;
            T(t,t)=sqrt(1-K(i,i));
            T(t,1:t-1)=(T(1:t-1,1:t-1)\K(B(1,1:b),i))';
        end
        b=b+1;
        B(1,b)=i;
    end
end
A=A(1,1:a);
end

