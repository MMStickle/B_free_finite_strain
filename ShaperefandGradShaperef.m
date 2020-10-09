function [Ninterp,GradNinterp]=ShaperefandGradShaperef(nnodes,ngauss)

     % funciones de interpolación y sus gradientes en el elemento de
     % referencia para elementos cuadriláteros de 4, 8 y 9 nodos,
     % con 4 y 9 puntos de integración. Ya evaluadas en los puntos de
     % integración.
     

   Ninterp=zeros(nnodes,ngauss);
   GradNinterp=zeros(2,nnodes,ngauss);
   
   if abs(nnodes-4)<1e-10  
     
     if abs(ngauss-4)<1e-10
         
            Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
                1/sqrt(3) -1/sqrt(3)
                1/sqrt(3) 1/sqrt(3)
                -1/sqrt(3) 1/sqrt(3)];
     
            for pto=1:ngauss   
  
                Ninterp(1,pto)=(1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))/4;
                Ninterp(2,pto)=(1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))/4;
                Ninterp(3,pto)=(1+Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))/4;
                Ninterp(4,pto)=(1-Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))/4;
                
                
                
                GradNinterp(1,1,pto)=(Gpd2con4(pto,2)-1)/4; GradNinterp(2,1,pto)=(Gpd2con4(pto,1)-1)/4;
                GradNinterp(1,2,pto)=(1-Gpd2con4(pto,2))/4; GradNinterp(2,2,pto)=-(Gpd2con4(pto,1)+1)/4;
                GradNinterp(1,3,pto)=(Gpd2con4(pto,2)+1)/4; GradNinterp(2,3,pto)=(Gpd2con4(pto,1)+1)/4;
                GradNinterp(1,4,pto)=-(Gpd2con4(pto,2)+1)/4; GradNinterp(2,4,pto)=(1-Gpd2con4(pto,1))/4;
                
                
          
            end    
        
     elseif abs(ngauss-9)<1e-10
         
         
            Gpd2con9=[0 0
                    -sqrt(3/5) -sqrt(3/5)
                        0 -sqrt(3/5)
                    sqrt(3/5) -sqrt(3/5)
                    sqrt(3/5) 0
                    sqrt(3/5) sqrt(3/5)
                        0 sqrt(3/5)
                    -sqrt(3/5) sqrt(3/5)
                    -sqrt(3/5) 0];
     
            for pto=1:ngauss   
  
                Ninterp(1,pto)=(1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))/4;
                Ninterp(2,pto)=(1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))/4;
                Ninterp(3,pto)=(1+Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))/4;
                Ninterp(4,pto)=(1-Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))/4;
                
                
                GradNinterp(1,1,pto)=(Gpd2con9(pto,2)-1)/4; GradNinterp(2,1,pto)=(Gpd2con9(pto,1)-1)/4;
                GradNinterp(1,2,pto)=(1-Gpd2con9(pto,2))/4; GradNinterp(2,2,pto)=-(Gpd2con9(pto,1)+1)/4;
                GradNinterp(1,3,pto)=(Gpd2con9(pto,2)+1)/4; GradNinterp(2,3,pto)=(Gpd2con9(pto,1)+1)/4;
                GradNinterp(1,4,pto)=-(Gpd2con9(pto,2)+1)/4; GradNinterp(2,4,pto)=(1-Gpd2con9(pto,1))/4;
                
                
                
       
            end    
        
     end
     
   elseif abs(nnodes-8)<1e-10
       
        if abs(ngauss-4)<1e-10
         
            Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
                1/sqrt(3) -1/sqrt(3)
                1/sqrt(3) 1/sqrt(3)
                -1/sqrt(3) 1/sqrt(3)];
     
            for pto=1:ngauss   
                
                Ninterp(1,pto)=(-(1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))*(1+Gpd2con4(pto,1)+Gpd2con4(pto,2))/4);
                Ninterp(2,pto)=(-(1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)+Gpd2con4(pto,2))/4);
                Ninterp(3,pto)=(-(1+Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)-Gpd2con4(pto,2))/4);
                Ninterp(4,pto)=(-(1-Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))*(1+Gpd2con4(pto,1)-Gpd2con4(pto,2))/4);
                Ninterp(5,pto)=((1-Gpd2con4(pto,1)^2)*(1-Gpd2con4(pto,2))/2);
                Ninterp(6,pto)=((1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2)/2);
                Ninterp(7,pto)=((1-Gpd2con4(pto,1)^2)*(1+Gpd2con4(pto,2))/2);
                Ninterp(8,pto)=((1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2)/2);
                
                
                
                
                
                GradNinterp(1,1,pto)=(1-Gpd2con4(pto,2))*(2*Gpd2con4(pto,1)+Gpd2con4(pto,2))/4; GradNinterp(2,1,pto)=(1-Gpd2con4(pto,1))*(2*Gpd2con4(pto,2)+Gpd2con4(pto,1))/4;
                GradNinterp(1,2,pto)=(Gpd2con4(pto,2)-1)*(Gpd2con4(pto,2)-2*Gpd2con4(pto,1))/4; GradNinterp(2,2,pto)=(1+Gpd2con4(pto,1))*(2*Gpd2con4(pto,2)-Gpd2con4(pto,1))/4;
                GradNinterp(1,3,pto)=(1+Gpd2con4(pto,2))*(2*Gpd2con4(pto,1)+Gpd2con4(pto,2))/4; GradNinterp(2,3,pto)=(1+Gpd2con4(pto,1))*(2*Gpd2con4(pto,2)+Gpd2con4(pto,1))/4;
                GradNinterp(1,4,pto)=(1+Gpd2con4(pto,2))*(2*Gpd2con4(pto,1)-Gpd2con4(pto,2))/4; GradNinterp(2,4,pto)=(Gpd2con4(pto,1)-1)*(Gpd2con4(pto,1)-2*Gpd2con4(pto,2))/4;
                GradNinterp(1,5,pto)=Gpd2con4(pto,1)*(Gpd2con4(pto,2)-1); GradNinterp(2,5,pto)=0.5*(Gpd2con4(pto,1)^2-1);
                GradNinterp(1,6,pto)=0.5*(1-Gpd2con4(pto,2)^2); GradNinterp(2,6,pto)=-Gpd2con4(pto,2)*(Gpd2con4(pto,1)+1);
                GradNinterp(1,7,pto)=-Gpd2con4(pto,1)*(Gpd2con4(pto,2)+1); GradNinterp(2,7,pto)=0.5*(1-Gpd2con4(pto,1)^2);
                GradNinterp(1,8,pto)=0.5*(Gpd2con4(pto,2)^2-1); GradNinterp(2,8,pto)=Gpd2con4(pto,2)*(Gpd2con4(pto,1)-1);
                
                         
                
  
              
            end    
        
        elseif abs(ngauss-9)<1e-10
         
         
            Gpd2con9=[0 0
                    -sqrt(3/5) -sqrt(3/5)
                        0 -sqrt(3/5)
                    sqrt(3/5) -sqrt(3/5)
                    sqrt(3/5) 0
                    sqrt(3/5) sqrt(3/5)
                        0 sqrt(3/5)
                    -sqrt(3/5) sqrt(3/5)
                    -sqrt(3/5) 0];
     
            for pto=1:ngauss   
                
                Ninterp(1,pto)=(-(1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))*(1+Gpd2con9(pto,1)+Gpd2con9(pto,2))/4);
                Ninterp(2,pto)=(-(1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)+Gpd2con9(pto,2))/4);
                Ninterp(3,pto)=(-(1+Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)-Gpd2con9(pto,2))/4);
                Ninterp(4,pto)=(-(1-Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))*(1+Gpd2con9(pto,1)-Gpd2con9(pto,2))/4);
                Ninterp(5,pto)=((1-Gpd2con9(pto,1)^2)*(1-Gpd2con9(pto,2))/2);
                Ninterp(6,pto)=((1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2)/2);
                Ninterp(7,pto)=((1-Gpd2con9(pto,1)^2)*(1+Gpd2con9(pto,2))/2);
                Ninterp(8,pto)=((1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2)/2);
                
                
                GradNinterp(1,1,pto)=(1-Gpd2con9(pto,2))*(2*Gpd2con9(pto,1)+Gpd2con9(pto,2))/4; GradNinterp(2,1,pto)=(1-Gpd2con9(pto,1))*(2*Gpd2con9(pto,2)+Gpd2con9(pto,1))/4;
                GradNinterp(1,2,pto)=(Gpd2con9(pto,2)-1)*(Gpd2con9(pto,2)-2*Gpd2con9(pto,1))/4; GradNinterp(2,2,pto)=(1+Gpd2con9(pto,1))*(2*Gpd2con9(pto,2)-Gpd2con9(pto,1))/4;
                GradNinterp(1,3,pto)=(1+Gpd2con9(pto,2))*(2*Gpd2con9(pto,1)+Gpd2con9(pto,2))/4; GradNinterp(2,3,pto)=(1+Gpd2con9(pto,1))*(2*Gpd2con9(pto,2)+Gpd2con9(pto,1))/4;
                GradNinterp(1,4,pto)=(1+Gpd2con9(pto,2))*(2*Gpd2con9(pto,1)-Gpd2con9(pto,2))/4; GradNinterp(2,4,pto)=(Gpd2con9(pto,1)-1)*(Gpd2con9(pto,1)-2*Gpd2con9(pto,2))/4;
                GradNinterp(1,5,pto)=Gpd2con9(pto,1)*(Gpd2con9(pto,2)-1); GradNinterp(2,5,pto)=0.5*(Gpd2con9(pto,1)^2-1);
                GradNinterp(1,6,pto)=0.5*(1-Gpd2con9(pto,2)^2); GradNinterp(2,6,pto)=-Gpd2con9(pto,2)*(Gpd2con9(pto,1)+1);
                GradNinterp(1,7,pto)=-Gpd2con9(pto,1)*(Gpd2con9(pto,2)+1); GradNinterp(2,7,pto)=0.5*(1-Gpd2con9(pto,1)^2);
                GradNinterp(1,8,pto)=0.5*(Gpd2con9(pto,2)^2-1); GradNinterp(2,8,pto)=Gpd2con9(pto,2)*(Gpd2con9(pto,1)-1);
                
                
                
                
                
                
                
       
            end    
        
        end
        
   elseif abs(nnodes-9)<1e-10
       
        if abs(ngauss-4)<1e-10
         
            Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
                1/sqrt(3) -1/sqrt(3)
                1/sqrt(3) 1/sqrt(3)
                -1/sqrt(3) 1/sqrt(3)];
     
            for pto=1:ngauss   
                
                Ninterp(1,pto)=(1/4)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2));
                Ninterp(2,pto)=(1/4)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2));
                Ninterp(3,pto)=(1/4)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2));
                Ninterp(4,pto)=(1/4)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2));
                
                Ninterp(5,pto)=(1/2)*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)^2);
                Ninterp(6,pto)=(1/2)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2);
                Ninterp(7,pto)=(1/2)*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)^2);
                Ninterp(8,pto)=(1/2)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2);
                              
                
                Ninterp(9,pto)=((1-Gpd2con4(pto,1)^2)*(1-Gpd2con4(pto,2)^2));
                
                
                
                
                
                GradNinterp(1,1,pto)=-(Gpd2con4(pto,1)/2 - 1/4)*(- Gpd2con4(pto,2)^2 + Gpd2con4(pto,2)); GradNinterp(2,1,pto)=-(2*Gpd2con4(pto,2) - 1)*(- Gpd2con4(pto,1)^2/4 + Gpd2con4(pto,1)/4);
                GradNinterp(1,2,pto)=-(Gpd2con4(pto,1)/2 + 1/4)*(- Gpd2con4(pto,2)^2 + Gpd2con4(pto,2)); GradNinterp(2,2,pto)=(2*Gpd2con4(pto,2) - 1)*(Gpd2con4(pto,1)^2/4 + Gpd2con4(pto,1)/4);
                GradNinterp(1,3,pto)=(Gpd2con4(pto,1)/2 + 1/4)*(Gpd2con4(pto,2)^2 + Gpd2con4(pto,2)); GradNinterp(2,3,pto)=(2*Gpd2con4(pto,2) + 1)*(Gpd2con4(pto,1)^2/4 + Gpd2con4(pto,1)/4);
                GradNinterp(1,4,pto)=(Gpd2con4(pto,1)/2 - 1/4)*(Gpd2con4(pto,2)^2 + Gpd2con4(pto,2)); GradNinterp(2,4,pto)=-(2*Gpd2con4(pto,2) + 1)*(- Gpd2con4(pto,1)^2/4 + Gpd2con4(pto,1)/4);
                GradNinterp(1,5,pto)=2*Gpd2con4(pto,1)*(- Gpd2con4(pto,2)^2/2 + Gpd2con4(pto,2)/2); GradNinterp(2,5,pto)=-(Gpd2con4(pto,1)^2 - 1)*(Gpd2con4(pto,2) - 1/2);
                GradNinterp(1,6,pto)=-(Gpd2con4(pto,2)^2 - 1)*(Gpd2con4(pto,1) + 1/2); GradNinterp(2,6,pto)=-2*Gpd2con4(pto,2)*(Gpd2con4(pto,1)^2/2 + Gpd2con4(pto,1)/2);
                GradNinterp(1,7,pto)=-2*Gpd2con4(pto,1)*(Gpd2con4(pto,2)^2/2 + Gpd2con4(pto,2)/2); GradNinterp(2,7,pto)=-(Gpd2con4(pto,1)^2 - 1)*(Gpd2con4(pto,2) + 1/2);
                GradNinterp(1,8,pto)=-(Gpd2con4(pto,2)^2 - 1)*(Gpd2con4(pto,1) - 1/2); GradNinterp(2,8,pto)=2*Gpd2con4(pto,2)*(- Gpd2con4(pto,1)^2/2 + Gpd2con4(pto,1)/2);
                GradNinterp(1,9,pto)=2*Gpd2con4(pto,1)*(Gpd2con4(pto,2)^2 - 1); GradNinterp(2,9,pto)=2*Gpd2con4(pto,2)*(Gpd2con4(pto,1)^2 - 1);
                         
                
  
              
            end    
        
        elseif abs(ngauss-9)<1e-10
         
         
            Gpd2con9=[0 0
                    -sqrt(3/5) -sqrt(3/5)
                        0 -sqrt(3/5)
                    sqrt(3/5) -sqrt(3/5)
                    sqrt(3/5) 0
                    sqrt(3/5) sqrt(3/5)
                        0 sqrt(3/5)
                    -sqrt(3/5) sqrt(3/5)
                    -sqrt(3/5) 0];
     
            for pto=1:ngauss   
  
                Ninterp(1,pto)=(1/4)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2));
                Ninterp(2,pto)=(1/4)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2));
                Ninterp(3,pto)=(1/4)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2));
                Ninterp(4,pto)=(1/4)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2));
                
                Ninterp(5,pto)=(1/2)*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)^2);
                Ninterp(6,pto)=(1/2)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2);
                Ninterp(7,pto)=(1/2)*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)^2);
                Ninterp(8,pto)=(1/2)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2);
                              
                
                Ninterp(9,pto)=((1-Gpd2con9(pto,1)^2)*(1-Gpd2con9(pto,2)^2));
                
                
                
                
                
                GradNinterp(1,1,pto)=-(Gpd2con9(pto,1)/2 - 1/4)*(- Gpd2con9(pto,2)^2 + Gpd2con9(pto,2)); GradNinterp(2,1,pto)=-(2*Gpd2con9(pto,2) - 1)*(- Gpd2con9(pto,1)^2/4 + Gpd2con9(pto,1)/4);
                GradNinterp(1,2,pto)=-(Gpd2con9(pto,1)/2 + 1/4)*(- Gpd2con9(pto,2)^2 + Gpd2con9(pto,2)); GradNinterp(2,2,pto)=(2*Gpd2con9(pto,2) - 1)*(Gpd2con9(pto,1)^2/4 + Gpd2con9(pto,1)/4);
                GradNinterp(1,3,pto)=(Gpd2con9(pto,1)/2 + 1/4)*(Gpd2con9(pto,2)^2 + Gpd2con9(pto,2)); GradNinterp(2,3,pto)=(2*Gpd2con9(pto,2) + 1)*(Gpd2con9(pto,1)^2/4 + Gpd2con9(pto,1)/4);
                GradNinterp(1,4,pto)=(Gpd2con9(pto,1)/2 - 1/4)*(Gpd2con9(pto,2)^2 + Gpd2con9(pto,2)); GradNinterp(2,4,pto)=-(2*Gpd2con9(pto,2) + 1)*(- Gpd2con9(pto,1)^2/4 + Gpd2con9(pto,1)/4);
                GradNinterp(1,5,pto)=2*Gpd2con9(pto,1)*(- Gpd2con9(pto,2)^2/2 + Gpd2con9(pto,2)/2); GradNinterp(2,5,pto)=-(Gpd2con9(pto,1)^2 - 1)*(Gpd2con9(pto,2) - 1/2);
                GradNinterp(1,6,pto)=-(Gpd2con9(pto,2)^2 - 1)*(Gpd2con9(pto,1) + 1/2); GradNinterp(2,6,pto)=-2*Gpd2con9(pto,2)*(Gpd2con9(pto,1)^2/2 + Gpd2con9(pto,1)/2);
                GradNinterp(1,7,pto)=-2*Gpd2con9(pto,1)*(Gpd2con9(pto,2)^2/2 + Gpd2con9(pto,2)/2); GradNinterp(2,7,pto)=-(Gpd2con9(pto,1)^2 - 1)*(Gpd2con9(pto,2) + 1/2);
                GradNinterp(1,8,pto)=-(Gpd2con9(pto,2)^2 - 1)*(Gpd2con9(pto,1) - 1/2); GradNinterp(2,8,pto)=2*Gpd2con9(pto,2)*(- Gpd2con9(pto,1)^2/2 + Gpd2con9(pto,1)/2);
                GradNinterp(1,9,pto)=2*Gpd2con9(pto,1)*(Gpd2con9(pto,2)^2 - 1); GradNinterp(2,9,pto)=2*Gpd2con9(pto,2)*(Gpd2con9(pto,1)^2 - 1);
                
                
                
                
                
                
                
       
            end    
        
        end
       
   end