function [Y] = spect(lambda,vect)
%   funkcja zwraca losowy zbiór z L-procesu przy danej dekompozycji
%   macierzy L
%   lamda- pionowy wektor wartości własnych
%   vect- macierz której kolumnami są kolejne wektory własne odpowaidającym
%   kolejnym lambdą

n=size(lambda,1);
if n~=size(vect,1) || n~=size(vect,2)
    disp('zły rozmiar danych');
    return ;
end
for i=1:n
    if abs(transpose(vect(:,i))*vect(:,i)-1)>1e-6
        disp('To nie jest baza ortogonalna');
        return;
    end
    for j=i+1:n
       if abs(transpose(vect(:,j))*vect(:,i))>1e-6
            disp('To nie jest baza ortonormalna');
            return;
       end
    end
end

J=[];
Y=[];
ile=0;   
e=eye(n);
for i=1:n %pierwsza pętla// wybieranie "aktywnych" wektorów własnych
    p=rand;
    if lambda(i)/(lambda(i)+1)>p
        J=union(J,i);
        ile=ile+1;
    end
end 
v=vect;
while ile>0 %druga pętla // wybieranie punktów do zbioru wynikowgo
    p=rand;%losowanie
    z=0;
    j=1;
    while z/ile<p %sprawdzanie który punkt wylosowaliśmy
        z=z+(v(j,J)*transpose(v(j,J)));
        j=j+1;
        if j==n+1 && z/ile<p
            r=0;
            for i=1:n
                r=r+v(i,J)*transpose(v(i,J));
            end
            p=p*sqrt(r);
            j=1;
        end
    end
    j=j-1;
    rzut=zeros(n,1);

    for x=J %ortogonalizacja bazy podprzestrzeni
        rzut=rzut+(e(j,:)*v(:,x))*v(:,x);  
    end
    rzut=rzut/sqrt(transpose(rzut)*rzut);
    for x=J
        v(:,x)=v(:,x)-((transpose(v(:,x))*rzut)*rzut);
        for y=J
            if y<x
            v(:,x)=v(:,x)-(transpose(v(:,x))*v(:,y))*v(:,y);
            end
        end
        if transpose(v(:,x))*v(:,x)>=1e-8
            v(:,x)=v(:,x)/sqrt(transpose(v(:,x))*v(:,x));
        else
             J=setdiff(J,x);
        end 
    end
    v(j,:)=zeros(1,n);
    Y=union(Y,j);
    ile=ile-1;
end

end
