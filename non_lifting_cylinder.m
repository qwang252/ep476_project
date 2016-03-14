%non-lifting flow around a cylinder of radius R with freestream
function non_lifting_cylinder(R,npanels,v_infinity)
v_freestream=v_infinity*[1 0 0];
%define the panels as well as end points using complex 
%Npanels represents the number of panels, and therefore defines the
%discricized shape of cylinder
dangle=2*pi/npanels;
for m=1:npanels
    theta(:,m)=+dangle*(m-1)-dangle/2;
end
%index of end points (xm,ym)
for m=1:npanels
    x(:,m)=real(R*exp(1i*theta(1,m)));
    y(:,m)=imag(R*exp(1i*theta(1,m)));
end
%define normal and tangetial directions for each panel
for m=1:npanels
    if m<=npanels-1
    t_hat(:,m)=[x(1,m+1)-x(1,m) y(1,m+1)-y(1,m) 0]/sqrt((x(1,m+1)-x(1,m))^2+(y(1,m+1)-y(1,m))^2);
    else
    t_hat(:,m)=[x(1,1)-x(1,m) y(1,1)-y(1,m) 0]/sqrt((x(1,1)-x(1,m))^2+(y(1,1)-y(1,m))^2); 
    end  
end
%good
z_hat=[0 0 1];
 for m=1:npanels
   n_hat(:,m)=cross(z_hat,t_hat(:,m));   
 end
%define midpoints for each panel
 for m=1:npanels
    if m<=npanels-1 
    midpt(:,m)=[x(1,m+1)+x(1,m) y(1,m+1)+y(1,m) 0]/2;
    else
    midpt(:,m)=[x(1,1)+x(1,m) y(1,1)+y(1,m) 0]/2;  
    end
 end 
%governing velocity equation
A=zeros(npanels);
T_velo=zeros(npanels);
 for i=1:npanels
     for j=1:npanels
        if i~=j
        if j<npanels 
        A(j,i)=(flat_panel_velocity(x(1,j),y(1,j),x(1,j+1),y(1,j+1),midpt(1,i),midpt(2,i),1))*n_hat(:,i);
        else
        A(j,i)=(flat_panel_velocity(x(1,j),y(1,j),x(1,1),y(1,1),midpt(1,i),midpt(2,i),1))*n_hat(:,i);
        end
        end      
     end
 end
lambda=zeros(1,npanels);
B=-v_freestream*n_hat;
lambda=B*inv(A-eye(npanels)/2)
%tangetial velocity
for i=1:npanels
     for j=1:npanels
        if i~=j
        if j<npanels 
        T_velo(j,i)=(flat_panel_velocity(x(1,j),y(1,j),x(1,j+1),y(1,j+1),midpt(1,i),midpt(2,i),lambda(1,j)))*t_hat(:,i);
        else
        T_velo(j,i)=(flat_panel_velocity(x(1,j),y(1,j),x(1,1),y(1,1),midpt(1,i),midpt(2,i),lambda(1,j)))*t_hat(:,i);
        end
        end      
     end
end
%tangetial velocity for each panel
T_velo
t_velo=sum(T_velo)+v_freestream*t_hat
%pressure coefficient 
for j=1:npanels+1
  if j<=npanels   
 c_p(:,j)=1-(t_velo(:,j)^2/v_infinity^2)
  else 
  c_p(:,j)=1-(t_velo(:,1)^2/v_infinity^2)
 end
end
%plot cp 
for i=1:npanels+1
    t(:,i)=dangle*(i-1);
end
tt=[0:1/2^6:1]*2*pi;
c_p_an=1-4*(sin(tt)).^2;
figure (1)
scatter(t,c_p)
legend('numerical solution');
 hold on,
plot(tt,c_p_an)
xlabel('theta','fontsize',16);
ylabel('Cp','fontsize',16);
%plot tangential velocity
v_tan_anly=2*v_infinity*sin(tt);
figure (2)
for i=1:npanels+1
    if i<=npanels
        t_velo_p(:,i)=t_velo(:,i);
    else
        t_velo_p(:,i)=t_velo(:,1);
    end

end
scatter(t,-t_velo_p),legend('numerical solution');
hold on
plot(tt,v_tan_anly);
end

 





 
 
 
 
 
 
 
 
 
 
 
 
 

