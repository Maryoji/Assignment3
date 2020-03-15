%%
% *FINITE DIFFERENCE METHOD*
%% QUESTION 2
% Using the Finite Difference Method in Assignment-2 to calculate the electric field and providing a field for the Monte-Carlo bottle-neck simulation.


W = 2;
L = 3;
V0 = 1;

dx = 0.025; % x mesh spacing
dy = 0.025; % y mesh spacing
nx = L/dx; % Number of points along x
ny = W/dy; % Number of points along y

Lb = 20;
Wb = 10;
% Generating the map of conductivity of the area
sigma_conduct = 1;
sigma_insulate = 10e-2;
% Construct the C matrix:
C = sigma_conduct.*ones(ny,nx);
Csubtract = zeros(ny,nx);

for x=1:nx
    for y=1:ny
        xx = x*dx;
        yy = y*dy;
        
        % The resistivity is made high in the rectangular regions:
        if(xx <= (L+Lb)/2 && xx >= (L-Lb)/2 && (yy >= W-Wb || yy <= Wb))
            Csubtract(y,x) = sigma_conduct-sigma_insulate;
        end
    end
end

% Filter the condicivity to avoid numerical issues that can occur if the derivatives are large.
Csubtract = imgaussfilt(Csubtract, 1);
C = C - Csubtract;

%%
% Below, the conducitivity is plotted. I performed some filtering/smoothing
% so that the derivative is not very large (approaching infinity) at the
% junction of the two regions. 

% figure(1);
% surf(linspace(0,L,nx),linspace(0,W,ny),C);
% title('Conductivity');
% view(30,45);
% xlabel('x (m)');
% ylabel('y (m)');
% grid on;

G = zeros(nx*ny,nx*ny);
F = zeros(nx*ny,1);

dx2 = 1./(dx.^2);
dy2 = 1./(dy.^2);

for x=2:(nx-1)
    for y=2:(ny-1)
        index = coordinate(x,y,nx);
        
        % Apply the equation derived earlier:
        G(index,index) = -2.*C(y,x).*(dx2 + dy2);
        G(index, coordinate(x+1,y,nx)) = dx2.*(0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        G(index, coordinate(x-1,y,nx)) = dx2.*(-0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        
        G(index, coordinate(x,y+1,nx)) = dy2.*(0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
        G(index, coordinate(x,y-1,nx)) = dy2.*(-0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
    end
end

%%
% Next, the F matrix is generated.

% The top and bottom boundaries
for x=2:(nx-1)
    index = coordinate(x,1,nx);
    G(index,index) = 1;
    G(index,coordinate(x,2,nx)) = -1;
    F(index) = 0;
    
    index = coordinate(x,ny,nx);
    G(index,index) = 1;
    G(index,coordinate(x,ny-1,nx)) = -1;
    F(index) = 0;
end

% The vertical boundaries
for y=1:ny
    index = coordinate(1,y,nx);
    G(index,index) = 1;
    F(index) = V0;
    
    index = coordinate(nx,y,nx);
    G(index,index) = 1;
    F(index) = 0;
end

%%
%
V = G\F;
V = reshape(V,[],ny)';

% figure(2);
% surf(linspace(0,L,nx),linspace(0,W,ny),V);
% view(30,45);
% xlabel('x (m)');
% ylabel('y (m)');
% title('Electric Potential (V)');
% grid on;

%%
% The electric field is $E=-\nabla V$. Here it is plotted along with the
% voltage.

%figure(3);
[Ex,Ey] = gradient(V,dx,dy);
Ex = -1.*Ex;
Ey = -1.*Ey;
% quiver(linspace(0,L,nx),linspace(0,W,ny),Ex,Ey,4);
% xlabel('x (m)');
% ylabel('y (m)');
% title('Electric Field (V/m)');
% axis([0 L 0 W]);
% grid on;


