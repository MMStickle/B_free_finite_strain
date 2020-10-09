function [bor1,bor2]=stimaNeuu(coord,t)


g31=0;
%g32=1e3*(cos(2*pi*t)-1);
g32=0;


Gpd1con3=[-sqrt(0.6);0;sqrt(0.6)];%ojo así puesto se trata del elemento -1 1.
Wd1con3=[5/9;8/9;5/9];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
z=zeros(6,1);
cp=zeros(2,1);%esto nos permite distinguir entre grado de libertad x e y
                  %en la linea siguiente deberíamos multiplicar por el
                  %jacobino evaluado en los puntos de integración pero
                  %como aparece tambien dividido en invtrasF, no ponemos
                  %ninguno, es decir, los cancelamos.


    
   
 for pto1=1:3                                      
  cp(1)=g31*norm([Gpd1con3(pto1)*(coord(2,1)+coord(2,2)-2*coord(2,3))-0.5*coord(2,1)+0.5*coord(2,2),Gpd1con3(pto1)*(coord(1,1)+coord(1,2)-2*coord(1,3))-0.5*coord(1,1)+0.5*coord(1,2)])*Wd1con3(pto1);            
  cp(2)=g32*norm([Gpd1con3(pto1)*(coord(2,1)+coord(2,2)-2*coord(2,3))-0.5*coord(2,1)+0.5*coord(2,2),Gpd1con3(pto1)*(coord(1,1)+coord(1,2)-2*coord(1,3))-0.5*coord(1,1)+0.5*coord(1,2)])*Wd1con3(pto1);
   
   
  z([1;3;5],1)=z([1;3;5],1)+[-(Gpd1con3(pto1)/2)*(1-Gpd1con3(pto1));(Gpd1con3(pto1)/2)*(1+Gpd1con3(pto1));1-Gpd1con3(pto1)^2]*cp(1);                     
  z([2;4;6],1)=z([2;4;6],1)+[-(Gpd1con3(pto1)/2)*(1-Gpd1con3(pto1));(Gpd1con3(pto1)/2)*(1+Gpd1con3(pto1));1-Gpd1con3(pto1)^2]*cp(2); 
  
 end
 

bor1=z([1;3;5],1);
bor2=z([2;4;6],1);

                       