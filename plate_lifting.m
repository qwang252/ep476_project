function plate_lifting(alpha,G)
% plate_lifting lifting flow around flat plate by conformal mapping.
%  plate_lifting(alpha, G) plots flow around flat plate at angle of
%  attack alpha with clockwise circulation Gamma with G= Gamma/(4 pi RV) 
%  R = cylinder radius,  V = velocity at infinity
%
%  Includes functions cylinderflow and joukowski at end 
%
% adapted from plate_inclined_conformal.m script
% 
% Fabian Waleffe, UW Madison
% EMA 521, Aerodynamics, Fall 2015

%V=1;     % Velocity at Infinity
R=1;     % Cylinder radius, use R=1/4 for chord = 1
r=R*[1:1/2^4:6]'; t=[0:1/2^6:1]*2*pi; z=r*exp(1i*t); % r, t=theta, z
w=cylinderflow(z,R,G); % complex potential for lifting flow 
Z=joukowski(z, R,alpha); % map cylinder to plate at alpha
X=real(Z); Y=imag(Z);

% PLOTTING 
hold off
contour(X,Y,imag(w),40), hold on  % streamlines in Z plane (X,Y)

% Plot cylinder or flat plate
zcy=R*exp(1i*pi*[0:1/2^5:2])'; Zcy=joukowski(zcy,R,alpha); % cylinder, flat plate
plot(real(Zcy),imag(Zcy),'k')

xlabel('X/R','fontsize',16), ylabel('Y/R','fontsize',16)
axis image, axis([-5*R 5*R -5*R 5*R])

% Plot stagnation streamlines
plotstaglines(r,R,alpha,G), commandwindow

return

%-------------------------------------------------------------------------
function w=cylinderflow(z,R,G) 
%Complex potential for lifting flow around cylinder, G=Gamma/(4 pi RV)
w=((z+R^2./z)+1i*2*G*R*log(z/R));
%-------------------------------------------------------------------------
function Z=joukowski(z,R,alpha) 
%Joukowski map: cylinder to flat plate at alpha clockwise
c2=R^2*exp(-2*1i*alpha);
Z=z+c2./z;
%-------------------------------------------------------------------------
 function  plotstaglines(r,R,alpha,G)
% Stagnation points: dw/dz=0
 if G <= 1
    sint=2*G*(r/R).*log(r/R)./(1-(r/R).^2);
    t=asin(sint); zs1=r.*exp(1i*t); zs2=r.*exp(1i*(pi-t));
    Zs1=joukowski(zs1,R,alpha);
    Zs2=joukowski(zs2,R,alpha);
    plot(real(Zs1),imag(Zs1),'k--')
    plot(real(Zs2),imag(Zs2),'k--')
end

 


