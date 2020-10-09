%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%Esto significa
%%%%%%%%%%%%%%% comentarios respecto a la programación determinantes.
clear all
clc
tic
global ax bx ay by TOL...
       rhos g...
       alphatemp1 alphatemp2 alphatemp3 alphatemp4 alphatemp5 alphatemp6 ...
       lambda mu
       
       

%load coor4n40elem.txt
%load conectivi4n40elem.txt
%load coor8n40elem.txt
%load conectivi8n40elem.txt
load coor9n40elem.txt
load conectivi9n40elem.txt


conectividad9=conectivi9n40elem(:,2:10)';
%filas de conectividad9= contenido numeración global/posición de 1 a 9 numeración local 
%columnas de conectividad9= elementos empleados en la discretización
coordenadas9n=coor9n40elem(:,[2,3])';
%filas de coordenadas9n= coordenas x e y de cada nodo (9 nodos por elemento)
%columnas de coordenadas9n= nodos de la discretización (numeración global) 9 nodos por elemento

%conectividad4=conectivi4n40elem(:,2:5)';
%filas de conectividad4= contenido numeración global/posición de 1 a 4 numeración local 
%columnas de conectividad4= elementos empleados en la discretización
%coordenadas4n=coor4n40elem(:,[2,3])';
%filas de coordenadas4n= coordenas x e y de cada nodo (4 nodos por elemento)
%columnas de coordenadas4n= nodos de la discretización (numeración global) 4 nodos por elemento

%conectividad8=conectivi8n40elem(:,2:9)';
%filas de conectividad8= contenido numeración global/posición de 1 a 8 numeración local 
%columnas de conectividad8= elementos empleados en la discretización
%coordenadas8n=coor8n40elem(:,[2,3])';
%filas de coordenadas8n= coordenas x e y de cada nodo (8 nodos por elemento)
%columnas de coordenadas8n= nodos de la discretización (numeración global) 8 nodos por elemento

nnodes=size(conectividad9,1);%número de nodos por elemento. En este caso 9 al emplear coordenadas9n
n9=size(coordenadas9n,2);%número de nodos cuando hay 9 nodos por elemento
pgausselem=9; % puntos de gauss por elemento
%n4=size(coordenadas4n,2);%número de nodos cuando hay 4 nodos por elemento
%n8=size(coordenadas8n,2);%número de nodos cuando hay 8 nodos por elemento
ne=size(conectividad9,2);%número de elementos en la malla
coordptog=ptogauss(coordenadas9n,conectividad9,ne,pgausselem);
nptgauss=length(coordptog(1,:));


%[Ninterpolant4nodes4ptgauss,GradNinterpolant4nodes4ptgauss]=ShaperefandGradShaperef(4,4);
%[Ninterpolant9nodes4ptgauss,GradNinterpolant9nodes4ptgauss]=ShaperefandGradShaperef(9,4);

%Funciones de forma y gradientes de las funciones de forma en el elemento de 
%referencia para elementos cuadrilateros con 9 nodos y 9 puntos de gauss por elemento
%Ninterp=zeros(nnodes,ngauss);GradNinterp=zeros(2,nnodes,ngauss);
[Ninterpolant9nodes9pgausselem,GradNinterpolant9nodes9pgausselem]=ShaperefandGradShaperef(nnodes,pgausselem);

%Fhat44=Fhatquad(coordenadas4n,conectividad4,ne,pgausselem,GradNinterpolant4nodes4ptgauss);
%Fhat94=Fhatquad(coordenadas9n,conectividad9,ne,pgausselem,GradNinterpolant9nodes4ptgauss);

% Calcula Fhat como elemento isoparamétrico X=Nhat_alpha*X_alpha a
% través de la expresión Fhat=tensorproduct(X_alpha,GradhatN_alpha), ya
% evaluado en cada punto de Gauss de cada elemento de la malla. El primer 9
% se refiere a elementos de 9 nodos el segundo 9 se refiere a 9 puntos de
% gauss.
Fhat99=Fhatquad(coordenadas9n,conectividad9,ne,nnodes,pgausselem,GradNinterpolant9nodes9pgausselem);





%geometría y malla
ax=min(coordenadas9n(1,:));bx=max(coordenadas9n(1,:));
ay=min(coordenadas9n(2,:));by=max(coordenadas9n(2,:));
tinic=0;tfin=10;ht=(tfin-tinic)/100;
TOL=1e-10;



%Condiciones de borde del esqueleto. 
%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%

Contouru1=find(abs(coordenadas9n(2,:)-ay)<1e-10);
Contouru2=find(abs(coordenadas9n(1,:)-bx)<1e-10);
Contouru4=find(abs(coordenadas9n(1,:)-ax)<1e-10);
Contouru3=find(abs(coordenadas9n(2,:)-by)<1e-10);

%Contourpw3=find(abs(coordenadas4n(2,:)-by)<1e-10);

%elementos que tienen un lado en el contorno Contouru3 y el contorno Contourpw3
%Esta información la podemos sacar directamente de GID
Contouruelem3=[39,40];%para la columna son los elementos 39 y 40. En este caso está puesto directamente
Contourpwelem3=Contouruelem3;

nodosenbordeloadu=zeros(3,size(Contouruelem3,2));
%nodosenbordeload nos dice de los nodos que forman el elemento que tiene un
%lado en el contorno cuales están efectivamente en el contorno. Y los
%ordena por: esquina_esquina_centro
nodosenbordeloadpw=zeros(2,size(Contourpwelem3,2));
%Esta información se podría sacar directamente de Gid

normales4=zeros(1,size(Contouruelem3,2));
%normales4 nos dice cual es la normal del elemento de referencia que se 
%transforma en el lado que está en el contorno 

for i=1:size(Contouruelem3,2)
nodosenbordeloadu(:,i)=conectividad9(find(abs(coordenadas9n(2,conectividad9(:,Contouruelem3(i)))-by)<1e-10),Contouruelem3(i));
%nodosenbordeloadpw(:,i)=conectividad4(find(abs(coordenadas4n(2,conectividad4(:,Contourpwelem3(i)))-by)<1e-10),Contourpwelem3(i));

if find(abs(coordenadas9n(2,conectividad9(:,Contouruelem3(i)))-by)<1e-10)==[1,2,5]
    normales4(i)=1;
elseif find(abs(coordenadas9n(2,conectividad9(:,Contouruelem3(i)))-by)<1e-10)==[2,3,6]
    normales4(i)=2;
elseif find(abs(coordenadas9n(2,conectividad9(:,Contouruelem3(i)))-by)<1e-10)==[3,4,7]
    normales4(i)=3;
elseif find(abs(coordenadas9n(2,conectividad9(:,Contouruelem3(i)))-by)<1e-10)==[1,4,8]
    normales4(i)=4;
  
end
end


normalespw4=normales4;

%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%las lineas que siguen son para visualizar malla, elementos de 9 nodos y pontos de Gauss.
%Este meshgrid es sólo para poder ver nodos y puntos de gauss. No tiene ninguna otra finalidad
%[X,Y]=meshgrid(ax:0.5:bx,ay:0.5:by);
%Z=zeros(size(X));
%mesh(X,Y,Z)
%axis equal;
%axis([ax-1 bx+1 ay-1 by+1]);
%aquí termina la creación de la malla para representar nodos y puntos de gauss encima de ella


%este primer bucle es para pasar de elemento a elemento 
%for i=1:ne
%i

%este segundo bucle es para representar los puntos de gauss
%for j=1:pgausselem
%hold on
%plot(coordptog(1,j+(i-1)*pgausselem),coordptog(2,j+(i-1)*pgausselem),'b.');
%pause
%end
%hold off

%este tercer bucle es para poder representar nodos de cada elemento
%for j=1:nnodes
%hold on
%plot(coordenadas9n(1,conectividad9(j,i)),coordenadas9n(2,conectividad9(j,i)),'r.');
%pause
%end
%hold off
%end

%estas representaciones no son para cada elemento si no para la geometría entera
%plot(coordenadas9n(1,setdiff(1:n9,unique(union(Contouru4,union(Contouru2,Contouru1))))),coordenadas9n(2,setdiff(1:n9,unique(union(Contouru4,union(Contouru2,Contouru1))))),'ro');
%plot(coordenadas9n(1,nodosenbordeloadu),coordenadas9n(2,nodosenbordeloadu),'r+');
%plot(coordenadas4n(1,nodosenbordeloadpw),coordenadas4n(2,nodosenbordeloadpw),'r+');
%plot(coordenadas4n(1,FreeNodespw),coordenadas4n(2,FreeNodespw),'b+');
%hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Array que permite el ensamblaje con dos grados de libertad por nodo
%de elemento cuadratico (n9)
IDu=zeros(2,n9);
IDu(1,:)=1+((1:n9)-1)*2;%ensamblaje de desplazamiento en x
IDu(2,:)=2+((1:n9)-1)*2;%ensamblaje de desplazamiento en y

%Array que permite ensamblaje con un grado de libertad por nodo de elemento
%bilineal (n4)
%IDp=zeros(1,n4);
%IDp(1,:)=1:n4;

%Inicializamos las matrices dGdu, dGdtheta, dHdu, dHdtheta
dGdu=zeros(2*n9,2*n9);
%dGdtheta=zeros(2*n9,n4);
%dHdu=zeros(n4,2*n9);
%dHdtheta=zeros(n4,n4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FreeDegfreeu=setdiff((1:2*n9)',unique(union(IDu(1,Contouru1)',union(IDu(2,Contouru1)',union(IDu(1,Contouru2)',IDu(1,Contouru4)','rows'),'rows'),'rows')));
FreeDegfreeu=setdiff((1:2*n9)',IDu(2,Contouru1)');

%FreeNodespw=setdiff((1:n4)',unique(Contourpw3)'); 


%valores de los coeficientes para la integración temporal
Beta=0.3025;ganma=0.6;
alphatemp1=(1/(Beta*ht^2));alphatemp2=(1/(Beta*ht));alphatemp3=((1-2*Beta)/(2*Beta));
alphatemp4=(ganma/(Beta*ht));alphatemp5=1-(ganma/Beta);alphatemp6=(1-(ganma/(2*Beta)))*ht;

%valores para los coeficientes del comportamiento constitutivo del esqueleto
lambda=1e6;
mu=5e5;  
%E=mu*(3*lambda+2*mu)/(lambda+mu)
%mu=lambda/(2*(lambda+mu))

%fracciónes de volumen iniciales y aceleracion de gravedad
%fis0=0.58;
%fiw0=0.42;
g=10;

%Densidades intrinsecas de referencia
%rhow=1000;%kg/m3
rhos=2700;%kg/m3
%rho=fis0*rhos+fiw0*rhow;

%Permeabilidad y compresibilidad del agua
%krlsat=1e-1;%m/s
%Kw=2.2e10;%N/m^2

%alphatuning para componente viscosa
%alphatuning=0;


%Abrimos y cerramos ficheros para guardar resultados numericos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Primary variables at nodes
Ufile = fopen('desp.asc','w');
Vfile=fopen('velocidades.asc','w');
Afile=fopen('aceleraciones.asc','w');
%Tfile=fopen('theta.asc','w');
%Tpuntofile=fopen('thetapunto.asc','w');
%Tdospuntosfile=fopen('thetadospuntos.asc','w');

fclose(Ufile);fclose(Vfile);fclose(Afile);
%fclose(Tfile);fclose(Tpuntofile);fclose(Tdospuntosfile);



%Secundary variables at gauss points
Piola1file=fopen('Piola1.asc','w');
fclose(Piola1file);
%Piola1efectfile=fopen('Piola1efect.asc','w');
%Piola1viscoefectfile=fopen('Piola1viscoefect.asc','w');
%fclose(Piola1efectfile);
%fclose(Piola1efectfile);


%Other variables
Deformationgradientfile = fopen('Deformationgradient.asc','w');
tiempofile=fopen('tiempo.asc','w');

fclose(Deformationgradientfile );

%Condiciones iniciale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Condiciones iniciales de las Primary and secondary variables 
U0=zeros(2*n9,1);
print_file('desp.asc', U0');
V0=zeros(2*n9,1);
print_file('velocidades.asc', V0');
A0=zeros(2*n9,1);
print_file('aceleraciones.asc', A0');

Piola10=zeros(3,3,nptgauss);
print_file('Piola1.asc', Piola10);

Deformationgradient0=zeros(3,3,nptgauss)+eye(3,3);
print_file('Deformationgradient.asc', Deformationgradient0);

print_file('tiempo.asc', tinic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variables auxiliares para actualizar valores con el tiempo
%las variables Z(:,1) ó que terminen en old se refieren al tiempo n, mientras que las variables
%Z(:,2) ó que terminen en new se refieren al tiempo n+1.

ZU=zeros(size(U0,1),2);increU=zeros(2*n9,2);
ZU(:,1)=U0;

ZV=zeros(size(V0,1),2);
ZV(:,1)=V0;

ZA=zeros(size(A0,1),2);
ZA(:,1)=A0;

ZPiola1old=Piola10;
ZPiola1new=ZPiola1old;

ZDeformgradold=Deformationgradient0;
ZDeformgradnew=ZDeformgradold;


ti=tinic;



while ti<tfin
    increFext=zeros(2*n9,1);
    

    for j=1:size(nodosenbordeloadu,2)
        %Es nodosenbordeload y no Contouru3 por que hay que resolver una
        %integral de borde que se ensambla. Flujo=Neuman=resolver integral
        %de borde
        [bt1,bt2]=stimaNeuu(coordenadas9n(:,nodosenbordeloadu(:,j)),ti+ht);
        increFext(IDu(1,nodosenbordeloadu(:,j)),1)=increFext(IDu(1,nodosenbordeloadu(:,j)),1)+bt1;
        increFext(IDu(2,nodosenbordeloadu(:,j)),1)=increFext(IDu(2,nodosenbordeloadu(:,j)),1)+bt2;
    end                                                   
              
    for ite=1:15
        
      
        G=-increFext(:,1);% con esta linea G=-increFext(:,1) estamos indicando  
                          % que el vector G se está evaluando en dos partes. 
                          % En la primera incorporamos las tensiones superficiales externas, 
                          % que no cambian dentro del bucle de N-R (dead load). La segunda parte se
                          % evaluan las integrales que si sufren modificación dentro
                          % del bucle de N-R. Esta segunda parte se evalua dentro
                          % del bucle de N-R. En la ecuación 4.62 de
                          % Wriggers (Nonlinear_finite_element) sería la
                          % segunda integral
                          
        
             
        
        for i=1:ne
        %este bucle para i desde 1 hasta ne sirve para:
        %1) calcular a nivel elemental a través de stimaG
        %2) ensamblar a través de 
        %    G(IDu(1:2,conectividad9(:,i)),1)=G(IDu(1:2,conectividad9(:,i)),1)+stimaG 
              
                                                                         
        G(IDu(1:2,conectividad9(:,i)),1)=G(IDu(1:2,conectividad9(:,i)),1)+...
            stimaG(Fhat99(:,:,(1:pgausselem)+(i-1)*pgausselem),Ninterpolant9nodes9pgausselem,GradNinterpolant9nodes9pgausselem,...
            ZPiola1new(:,:,(1:pgausselem)+(i-1)*pgausselem),ZU(IDu(1:2,conectividad9(:,i)),2),...
            ZU(IDu(1:2,conectividad9(:,i)),1),ZV(IDu(1:2,conectividad9(:,i)),1),ZA(IDu(1:2,conectividad9(:,i)),1),...
                pgausselem,nnodes,ti+ht);          
          
                      
                                    
        end
        
       
        RHS=G(FreeDegfreeu,1);
        
       
        if ite==1
              norm0=norm(RHS);
        end
        fprintf('%i\t%i\t%i\t%i\t%i\n',ti,ht,ite,norm(RHS),norm(RHS)/norm0);
        
               
        if (norm(RHS)<TOL*norm0 || norm(RHS)<1e-7) 
                break
        end
               
        dGdu=zeros(2*n9,2*n9);
    
        
        
        for i=1:ne
        %este bucle para i desde 1 hasta ne sirve para:
        %1) calcular la matriz tangente a nivel elemental a través de stimadGdu
        %2) ensamblar a través de 
        %  dGdu(IDu(1:2,conectividad9(:,i)),IDu(1:2,conectividad9(:,i)))=...
        %  dGdu(IDu(1:2,conectividad9(:,i)),IDu(1:2,conectividad9(:,i)))+stimadGdu
        
        dGdu(IDu(1:2,conectividad9(:,i)),IDu(1:2,conectividad9(:,i)))=dGdu(IDu(1:2,conectividad9(:,i)),IDu(1:2,conectividad9(:,i)))+...
            stimadGdu(Fhat99(:,:,(1:pgausselem)+(i-1)*pgausselem),Ninterpolant9nodes9pgausselem,GradNinterpolant9nodes9pgausselem,...
            ZDeformgradnew(:,:,(1:pgausselem)+(i-1)*pgausselem),pgausselem,nnodes);
    
        end
        
        
        Dr=dGdu(FreeDegfreeu,FreeDegfreeu);
        
        
           
        aux=Dr\(-RHS);
        
           
        ZU(FreeDegfreeu,2)=ZU(FreeDegfreeu,2)+aux(:,1);
        
        
        for i=1:ne
           
            
            [Piola1,Deformationgradient]=updatesecondaryvariables(Fhat99(:,:,(1:pgausselem)+(i-1)*pgausselem),GradNinterpolant9nodes9pgausselem,...
                                                                                   ZU(IDu(1:2,conectividad9(:,i)),2),pgausselem,nnodes);
            
            
            ZPiola1new(1:2,1:2,(1:pgausselem)+(i-1)*pgausselem)=Piola1;
            ZDeformgradnew(1:2,1:2,(1:pgausselem)+(i-1)*pgausselem)=Deformationgradient;
                     
            
        end
        
        
    end    
    
    if abs(ite-15)<1e-10
        error('max ite alcanzado');
       
    end
    
        ZA(FreeDegfreeu,2)=alphatemp1*(ZU(FreeDegfreeu,2)-ZU(FreeDegfreeu,1))-alphatemp2*ZV(FreeDegfreeu,1)-alphatemp3*ZA(FreeDegfreeu,1);

        ZV(FreeDegfreeu,2)=alphatemp4*(ZU(FreeDegfreeu,2)-ZU(FreeDegfreeu,1))+alphatemp5*ZV(FreeDegfreeu,1)+alphatemp6*ZA(FreeDegfreeu,1);
           
   
      ti=ti+ht;
      ZU(:,1)=ZU(:,2);
      ZV(:,1)=ZV(:,2);
      ZA(:,1)=ZA(:,2);
      ZPiola1old(:,:,:)=ZPiola1new(:,:,:); 
      ZDeformgradold(:,:,:)=ZDeformgradnew(:,:,:);
     
      
     
    
       
      print_file('tiempo.asc',ti);
      print_file('desp.asc', ZU(:,2)')
      print_file('velocidades.asc', ZV(:,2)')
      print_file('aceleraciones.asc', ZA(:,2)');
      print_file('Piola1.asc', ZPiola1new(:,:,:));
      print_file('Deformationgradient.asc', ZDeformgradnew(:,:,:));
      
  
   
   %OJO:Zsigefectnew(:,:,:) no se anula tras pasar la información a 
   %Zsigefectold(:,:,:) ya que estamos considerando la fuerza externa de forma acumulativa y el
   %siguiente incremento de fuerza externa parte del estado de equilibrio
   %encontrado al final del  incremento anterior de fuerza externa.
      
      
      
     
        
end
 
load desp.asc
load tiempo.asc

%tn=linspace(1,size(desp,1),2);
for tn=1:size(desp,1)

for i=1:size(conectividad9,2)
     ord=[1,5,2,6,3,7,4,8,1];
     for j=1:9
        xnode(tn,j+(i-1)*9)=coordenadas9n(1,conectividad9(ord(j),i))+desp(tn,IDu(1,conectividad9(ord(j),i)));
        ynode(tn,j+(i-1)*9)=coordenadas9n(2,conectividad9(ord(j),i))+desp(tn,IDu(2,conectividad9(ord(j),i)));
        
     end
        
                
    
end
        plot(xnode(tn,:),ynode(tn,:),'b')
        hold on
        plot(xnode(tn,10:18),ynode(tn,10:18),'r.');
        hold off
        axis equal;
        axis([ax-5 bx+5 ay-2 by+2]);
        
       pause(0.1)
end

 


 


%Para poder ver la evolución de las 6 componentes de tensión en el punto de gauss ptogresul
%ptogresul=1;
%aa=1+(ptogresul-1)*3;
%sigefect(:,aa:(aa+2));
%para poder ver las 6 comoponentes de tensiones en toda la malla para un
%tiempo fijo tn
%n=1305;
%aa=1+(n-1)*3;
%sigefect(aa:(aa+2),1:3);

%para el primer punto de gauss en el segundo instante de tiempo sigefect(4:6,1:3)

%para ver la tension vertical del primer punto de gauss a lo largo del
%tiempo
%plot(tiempo,sigefect(2+(0:(length(tiempo)-1))*3,2))

%para ver el desplazamiento vertical de uno de los nodos superiores:
plot(tiempo,desp(:,IDu(2,find(abs(coordenadas9n(1,:)-1)<1e-3 & abs(coordenadas9n(2,:)-10)<1e-3))),'r')









%for i=1:size(desp,1)
%        totaldisp=sqrt(desp(i,IDu(1,unique(conectividad9([1,2,3],:)))).^2+desp(i,IDu(2,unique(conectividad9([1,2,3],:)))).^2);
%        coorplot(1,:)=coordenadas4n(1,:)+desp(i,IDu(1,unique(conectividad9([1,2,3],:))));      
%        coorplot(2,:)=coordenadas4n(2,:)+desp(i,IDu(2,unique(conectividad9([1,2,3],:))));
%        trisurf(conectividad4',coorplot(1,:)',coorplot(2,:)',totaldisp,'facecolor','interp','LineStyle','-','facelighting','phong')
%        title('total displacement [m]','FontSize',20,'FontWeight','bold')
%        colorbar('location','EastOutside','FontSize',20)
%        axis([ax-1 bx+1 ay-1 by+1]);
%        axis equal
%        view(2);
%        pause
              
% end





