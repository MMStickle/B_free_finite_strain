function Fhat=Fhatquad(coordenadas,conectividad,ne,nnodes,ngauss,GradNinterpolant)
% Calcula Fhat como elemento isoparam�trico X=Nhat_alpha*X_alpha a
% trav�s de la expresi�n Fhat=tensorproduct(X_alpha,GradhatN_alpha), ya
% evaluado en cada punto de Gauss de cada elemento de la malla.
Fhat=zeros(2,2,ngauss*ne);  
 for i=1:ne
     
   coord=coordenadas(:,conectividad(:,i));
      
   for pto=1:ngauss   
                 
     for nn=1:nnodes
         Fhat(:,:,pto+(i-1)*ngauss)=Fhat(:,:,pto+(i-1)*ngauss)+kron(coord(:,nn),GradNinterpolant(:,nn,pto)');
     end
    
   end
 end
            
    