function [P,F]=updatesecondaryvariables(Fhat,GradNref,unmas1i,pgausselem,nnodes)

% la i que acompaña a unmas1 en las variables de entrada viene a
% recordar que los valores de entrada Piola1i y unmas1i, para esta función son los de la
% iteración i de N-R para el tiempo n+1..

global lambda mu 


    
    
    F=zeros(2,2,pgausselem);
    P=zeros(2,2,pgausselem);
    
    
    
       
   
    
    for pto=1:pgausselem 
    
           
     Fhpto=Fhat(:,:,pto);
     gradNptog=Fhpto'\GradNref(:,:,pto);
     F(1:2,1:2,pto)=eye(2);
     for nn=1:nnodes
                
         F(1:2,1:2,pto)=F(1:2,1:2,pto)+kron(unmas1i((1:2)+(nn-1)*2,1),gradNptog(:,nn)');
         
     end
     
     Jpto=det(F(1:2,1:2,pto));
    %comportamiento elastico lineal
    P(1:2,1:2,pto)=(lambda/2)*(Jpto^2-1)*inv(F(1:2,1:2,pto)')+mu*(F(1:2,1:2,pto)-inv(F(1:2,1:2,pto)'));
     
    
    
    end
    
    


