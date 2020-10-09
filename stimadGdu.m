function dGdu=stimadGdu(Fhat,Nref,GradNref,Deformationgradient,pgausselem,nnodes)

                          

global lambda mu rhos alphatemp1 

%Si empleamos 9 puntos de integración son estos los pesos asociados a cada
%punto de integración.
%Centro-esquina-medio-esquina-medio-esquina-medio-esquina-medio
Wd2con9=[64/81;25/81;40/81;25/81;40/81;25/81;40/81;25/81;40/81];   
  
   z1=zeros(nnodes*2,nnodes*2);z2=z1;    
     
    for pto=1:pgausselem 
     %Deformationgradient=zeros(3,3,nptgauss)
     F=Deformationgradient(1:2,1:2,pto);Jpto=det(F);
     %GradNref=zeros(2,nnodes,ngauss);
     Fhpto=Fhat(:,:,pto);Jhpto=det(Fhpto);
     gradNptog=Fhpto'\GradNref(:,:,pto);
 
     
     
     for nn1=1:nnodes
         for nn2=1:nnodes
             
           z1((1:2)+(nn1-1)*2,(1:2)+(nn2-1)*2)=z1((1:2)+(nn1-1)*2,(1:2)+(nn2-1)*2)+...
               (lambda*Jpto^2*kron(F'\gradNptog(:,nn1),(F'\gradNptog(:,nn2))')+...
               mu*dot(gradNptog(:,nn1),gradNptog(:,nn2))*eye(2)+...
               (mu-(lambda/2)*(Jpto^2-1))*kron(F'\gradNptog(:,nn2),(F'\gradNptog(:,nn1))'))*Jhpto*Wd2con9(pto);
           
           z2((1:2)+(nn1-1)*2,(1:2)+(nn2-1)*2)=z2((1:2)+(nn1-1)*2,(1:2)+(nn2-1)*2)+...
               alphatemp1*rhos*Nref(nn1,pto)*Nref(nn2,pto)*Jhpto*Wd2con9(pto)*eye(2);
         end
     end
     
 
    
     
     
     
     
    end
    
    
dGdu=z1+z2;

