function [Y] = wizu(lambda,vect)
%   funkcja przedstawia kt�re punkty zosta�y wybrane a kt�re nie w
%   L-procesie z macierz� L, kt�rej dekompozycja na wektory w�asne jest
%   dana przez lambda i vect
%   lamda- pionowy wektor warto�ci w�asnych
%   vect- macierz kt�rej kolumnami s� kolejne wektory w�asne 

n=size(lambda,1);
if n~=size(vect,1) || n~=size(vect,2)
    disp('z�y rozmiar danych');
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

for i=1:n %pierwsza p�tla
    p=rand;
    if lambda(i)/(lambda(i)+1)>p
        J=union(J,i);
        ile=ile+1;
    end
end 
v=vect;
while ile>0 %druga p�tla
    p=rand;%losowanie
    z=0;
    j=1;
    while z/ile<p %sprawdzanie kt�ry punkt wylosowali�my
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
s=floor(sqrt(n));
x=zeros(n,1);
y=zeros(n,1);
for i=1:n
    x(i)=mod(i-1,s);
    y(i)=(i-1-mod(i-1,s))/s;
end
c = zeros(n,3);
for j=Y
    c(j,:)=[1,0,0];
end
U=setdiff(1:n,Y);
hold on
scatter(x(Y),y(Y),120,c(Y,:),'filled')
scatter(x(U),y(U),100,c(U,:),'filled')
hold off
axis off
legend({'\fontsize{14}Wybrane punkty','\fontsize{14}Niewybrane punkty'},'Location','southwestoutside');
end
