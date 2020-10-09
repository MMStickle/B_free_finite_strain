function Gue=stimaG(Fhat,Nref,GradNref,Piola1i,unmas1i,un,vn,an,pgausselem,nnodes,timas1)
      
global rhos g  alphatemp1 alphatemp2 alphatemp3
% la i que acompaña a Piola1 y a unmas1 en las variables de entrada viene a
% recordar que los valores de entrada Piola1i y unmas1i, para esta función son los de la
% iteración i de N-R.

%Si empleamos 9 puntos de integración son estos los pesos asociados a cada
%punto de integración.
%Centro-esquina-medio-esquina-medio-esquina-medio-esquina-medio
    Wd2con9=[64/81;25/81;40/81;25/81;40/81;25/81;40/81;25/81;40/81];   
    
%El residuo G que depende del desplazamiento u a nivel elemental e (Gue)
%tiene tantas filas como número de nodos en cada elemento por grados de
%libertad del desplazamiento
    Gue=zeros(nnodes*2,1);
    
    for pto=1:pgausselem 
    %Nref=zeros(nnodes,ngauss)
    anmas1ptog=zeros(2,1);
        for nn=1:nnodes
            anmas1ptog=anmas1ptog+Nref(nn,pto)*(alphatemp1*(unmas1i((1:2)+(nn-1)*2,1)-un((1:2)+(nn-1)*2,1))-...
                       alphatemp2*vn((1:2)+(nn-1)*2,1)-alphatemp3*an((1:2)+(nn-1)*2,1));
        end
   
    
    %%%%%%%%%%%%%%%%%
    %Gradiente de las funciones de forma
    %GradNref=zeros(2,nnodes,ngauss);
     Fhpto=Fhat(:,:,pto);Jhpto=det(Fhpto);
     gradNptog=Fhpto'\GradNref(:,:,pto);
     
     for nn=1:nnodes
     
     Gue((1:2)+(nn-1)*2,1)=Gue((1:2)+(nn-1)*2,1)+(Piola1i(1:2,1:2,pto)*gradNptog(:,nn)-...
         rhos*Nref(nn,pto)*([0;-min(timas1,1)*g]-anmas1ptog))*Jhpto*Wd2con9(pto);
   
     end
    end
       
 


