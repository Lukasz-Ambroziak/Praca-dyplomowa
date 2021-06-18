function [Y] = smac(lambda,vect)
%   tworzy kolorowa macierz z prawdopodobieñstwami 
%   wybrania kolejnych elementów w algorytmie spectralnym
%   lamda- pionowy wektor wartoœci w³asnych
%   vect- macierz której kolumnami s¹ kolejne wektory w³asne 

n=size(lambda,1);
if n~=size(vect,1) || n~=size(vect,2)
    disp('z³y rozmiar danych');
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

for i=1:n %pierwsza pêtla
    p=rand;
    if lambda(i)/(lambda(i)+1)>p
        J=union(J,i);
        ile=ile+1;
    end
end 
v=vect;
o=ile;
M=zeros(ile,n);
while ile>0 %druga pêtla
    p=rand;%losowanie
    z=0;
    j=1;
    while z/ile<p %sprawdzanie który punkt wylosowaliœmy
        z=z+(v(j,J)*transpose(v(j,J)));
        j=j+1;
        if j==n+1 && z/ile<p
            for i=1:n
                r=v(i,J)*transpose(v(i,J));
            end
            p=p*sqrt(r);
            j=1;
        end
    end
    j=j-1;
    
    for q=1:n
        a(q)=(v(q,J)*transpose(v(q,J)))/ile;
    end
    M(o-ile+1,:)=a;
    
    
    rzut=zeros(n,1);
    e=eye(n);
    for x=J %ortogonalizacja bazy podprzestrzeni
        rzut=rzut+(e(j,:)*v(:,x))*v(:,x);  
    end
    rzut=rzut/sqrt(transpose(rzut)*rzut);
    for x=J
        v(:,x)=v(:,x)-((transpose(v(:,x))*rzut)*rzut);
        for y=intersect(J,(1:x-1))
            v(:,x)=v(:,x)-(transpose(v(:,x))*v(:,y))*v(:,y);
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
heatmap(M)
colormap parula
end
