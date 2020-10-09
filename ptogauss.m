function punto=ptogauss(coord,nodelement,ne,pgausselem)

%Empleamos elementos isoparamétricos para establecer la posición de los
%puntos de gauss X(ptog)=Nhat_alpha(ptogr)*X_alpha donde 
%X(ptog) son las coordenadas del punto de gauss en el elemento físico en la configuración de referencia
%Nhat_alpha(ptog) son las funciones de forma en el elemento de referencia(isoparamétrico) evaluadas en las posiciones de los puntos de gauss en ese elemento de referencia.
%X_alpha son las coordenadas de los nodos de cada elemento físico en la configuración de referencia

punto=zeros(2,pgausselem*ne);
%filas de punto= coordenas x e y de cada punto de gauss
%columnas de punto= puntos de gauss en toda la malla. Se van rellenando
%elemento a elemento. 
%Primer elmento -> todos los puntos de gauss
%Segundo elemento -> todos los puntos de gauss
%y así sucesivamente

if abs(size(nodelement,1)-6)<1e-10
%puntos de gauss para elementos de 6 nodos
    
if pgausselem==3 %tres puntos de gauss por elemento

    
 Gpd2con3=[1/6 1/6
          2/3 1/6
          1/6 2/3];
   
    
for i=1:ne
    for pto=1:pgausselem 
    % X(ptog)=Nhat_alpha(ptogr)*X_alpha    
    punto(1:2,pto+(i-1)*pgausselem)=(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))*(2*(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))-1)*coord(:,nodelement(1,i))+...
                                    Gpd2con3(pto,1)*(2*Gpd2con3(pto,1)-1)*coord(:,nodelement(2,i))+...
                                    Gpd2con3(pto,2)*(2*Gpd2con3(pto,2)-1)*coord(:,nodelement(3,i))+...
                                    4*Gpd2con3(pto,1)*(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))*coord(:,nodelement(4,i))+...
                                    4*Gpd2con3(pto,1)*Gpd2con3(pto,2)*coord(:,nodelement(5,i))+...
                                    4*Gpd2con3(pto,2)*(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))*coord(:,nodelement(6,i));
                               
         
    end  
end

elseif pgausselem==7

    a=(6+sqrt(15))/21;
    b=(4/7)-a;

    Gpd2con7=[1/3 1/3
         a a
         1-2*a a
         a 1-2*a
         b b
         1-2*b b
         b 1-2*b];
   for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))*(2*(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))-1)*coord(:,nodelement(1,i))+...
                                    Gpd2con7(pto,1)*(2*Gpd2con7(pto,1)-1)*coord(:,nodelement(2,i))+...
                                    Gpd2con7(pto,2)*(2*Gpd2con7(pto,2)-1)*coord(:,nodelement(3,i))+...
                                    4*Gpd2con7(pto,1)*(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))*coord(:,nodelement(4,i))+...
                                    4*Gpd2con7(pto,1)*Gpd2con7(pto,2)*coord(:,nodelement(5,i))+...
                                    4*Gpd2con7(pto,2)*(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))*coord(:,nodelement(6,i));
             
    end
   end
end

elseif abs(size(nodelement,1)-4)<1e-10

 if pgausselem==4
 Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) 1/sqrt(3)
          -1/sqrt(3) 1/sqrt(3)];
   
    
for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=((1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    ((1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    ((1+Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    ((1-Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))/4)*coord(:,nodelement(4,i));
                                    
                                
                            
         
    end  
end

elseif pgausselem==9

    

    Gpd2con9=[0 0
          -sqrt(3/5) -sqrt(3/5)
              0 -sqrt(3/5)
           sqrt(3/5) -sqrt(3/5)
           sqrt(3/5) 0
           sqrt(3/5) sqrt(3/5)
           0 sqrt(3/5)
           -sqrt(3/5) sqrt(3/5)
           -sqrt(3/5) 0];
       
   for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=((1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    ((1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    ((1+Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    ((1-Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))/4)*coord(:,nodelement(4,i));
             
    end
   end
 end 
 
 elseif abs(size(nodelement,1)-8)<1e-10

 if pgausselem==4
 Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) 1/sqrt(3)
          -1/sqrt(3) 1/sqrt(3)];
   
    
for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=(-(1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))*(1+Gpd2con4(pto,1)+Gpd2con4(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    (-(1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)+Gpd2con4(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    (-(1+Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)-Gpd2con4(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    (-(1-Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))*(1+Gpd2con4(pto,1)-Gpd2con4(pto,2))/4)*coord(:,nodelement(4,i))+...
                                    ((1-Gpd2con4(pto,1)^2)*(1-Gpd2con4(pto,2))/2)*coord(:,nodelement(5,i))+...
                                    ((1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2)/2)*coord(:,nodelement(6,i))+...
                                    ((1-Gpd2con4(pto,1)^2)*(1+Gpd2con4(pto,2))/2)*coord(:,nodelement(7,i))+...
                                    ((1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2)/2)*coord(:,nodelement(8,i));
                               
         
    end  
end

elseif pgausselem==9

    

    Gpd2con9=[0 0
          -sqrt(3/5) -sqrt(3/5)
              0 -sqrt(3/5)
           sqrt(3/5) -sqrt(3/5)
           sqrt(3/5) 0
           sqrt(3/5) sqrt(3/5)
           0 sqrt(3/5)
           -sqrt(3/5) sqrt(3/5)
           -sqrt(3/5) 0];
       
   for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=(-(1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))*(1+Gpd2con9(pto,1)+Gpd2con9(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    (-(1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)+Gpd2con9(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    (-(1+Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)-Gpd2con9(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    (-(1-Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))*(1+Gpd2con9(pto,1)-Gpd2con9(pto,2))/4)*coord(:,nodelement(4,i))+...
                                    ((1-Gpd2con9(pto,1)^2)*(1-Gpd2con9(pto,2))/2)*coord(:,nodelement(5,i))+...
                                    ((1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2)/2)*coord(:,nodelement(6,i))+...
                                    ((1-Gpd2con9(pto,1)^2)*(1+Gpd2con9(pto,2))/2)*coord(:,nodelement(7,i))+...
                                    ((1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2)/2)*coord(:,nodelement(8,i));
             
    end
   end
 end   
elseif abs(size(nodelement,1)-9)<1e-10

 if pgausselem==4
 Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) 1/sqrt(3)
          -1/sqrt(3) 1/sqrt(3)];
   
    
for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=((1/4)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2)))*coord(:,nodelement(1,i))+...
                                    ((1/4)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2)))*coord(:,nodelement(2,i))+...
                                    ((1/4)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2)))*coord(:,nodelement(3,i))+...
                                    ((1/4)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2)))*coord(:,nodelement(4,i))+...
                                    ((1/2)*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)^2))*coord(:,nodelement(5,i))+...
                                    ((1/2)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2))*coord(:,nodelement(6,i))+...
                                    ((1/2)*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)^2))*coord(:,nodelement(7,i))+...
                                    ((1/2)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2))*coord(:,nodelement(8,i))+...
                                    (((1-Gpd2con4(pto,1)^2)*(1-Gpd2con4(pto,2)^2)))*coord(:,nodelement(9,i));
                               
                     
                                
                                
                                
    end  
end

elseif pgausselem==9

    

    Gpd2con9=[0 0
          -sqrt(3/5) -sqrt(3/5)
              0 -sqrt(3/5)
           sqrt(3/5) -sqrt(3/5)
           sqrt(3/5) 0
           sqrt(3/5) sqrt(3/5)
           0 sqrt(3/5)
           -sqrt(3/5) sqrt(3/5)
           -sqrt(3/5) 0];
       
   for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=((1/4)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2)))*coord(:,nodelement(1,i))+...                                     
                                    ((1/4)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2)))*coord(:,nodelement(2,i))+...                                     
                                    ((1/4)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2)))*coord(:,nodelement(3,i))+...                                     
                                    ((1/4)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2)))*coord(:,nodelement(4,i))+...                                   
                                    ((1/2)*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)^2))*coord(:,nodelement(5,i))+...                                     
                                    ((1/2)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2))*coord(:,nodelement(6,i))+...                                    
                                    ((1/2)*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)^2))*coord(:,nodelement(7,i))+...                                     
                                    ((1/2)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2))*coord(:,nodelement(8,i))+...                                     
                                    (((1-Gpd2con9(pto,1)^2)*(1-Gpd2con9(pto,2)^2)))*coord(:,nodelement(9,i));
                                  
                                
                                 
          

                
          
             
    end
   end
 end    
 
end

