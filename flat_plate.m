
%velocity at (x,y) induced by a flat source panel from (x1,y1) to (x2,y2) of
%strength lambda_1
function flat_plate(x1,y1,x2,y2,x,y,lambda_1)
%flat panel end point 1
p1=[x1 y1 0];
%flat panel end point 2
p2=[x2 y2 0];
%point of interested velocity 
p=[x y 0];
r1=p-p1;
r2=p-p2;
r1_mag=norm(r1);
r2_mag=norm(r2);
p1p2_mag=norm(p2-p1);
t_hat=1/p1p2_mag*(p2-p1)
x_hat=[1 0 0];
y_hat=[0 1 0];
z_hat=cross(x_hat,y_hat);
r1r2_cross=cross(r1,r2);
alpha=atan2(dot(r1r2_cross,z_hat),dot(r1,r2))
n_hat=cross(z_hat,t_hat);
flat_plate=lambda_1/2/pi*(log(r1_mag/r2_mag)*t_hat+alpha*n_hat)
% flat_plate



