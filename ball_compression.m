npoints = 20;

dx = 1/npoints;
radius = 1;
tfinal = 100;
nsteps = 10000;
dt = tfinal/nsteps;
g = 10;
tau = 1;
m = 1;
U = 0; %elastic potential
pos = 5;
vel = -5;

E_ball = 10^5; %Young's modulus of the ball's material (rubber)
E_wall = 23; %Young's modulus of the wall's material (concrete)
G_ball = .0006; %Shear modulus of the ball's material (rubber)
G_wall = .8; %Shear modulus of the wall's material

E_effective = 1/((1-G_ball^2)/E_ball + (1-G_wall^2)/E_wall);

delta_k = E_effective*dx; %stiffness of spring over area element

KE = 0;

[x,y,z] = create_ball(npoints, radius);
z = z + pos;

for t = 1:nsteps
    clf
    
    KE = .5*m*vel^2;
    normal_force = get_force(pos, E_effective, radius);
    a = normal_force/m - g;
    vel = vel + a*dt;
%   U = U + normal_force*dt;
%   KE = KE - normal_force*dt;
%   vel = sign*sqrt(2*abs(KE)/m) - g*dt;
    pos = pos + vel*dt;
    z = z + vel*dt;
        
    surfl(x,y,z); hold all
        axis([-5,5,-5,5,0,10])
    if mod(t, 2) == 0
    drawnow
    end
end    

function [Xcoords, Ycoords, Zcoords] = create_ball(npoints, radius)

    Xcoords=zeros(npoints, npoints);
    Ycoords=zeros(npoints, npoints);
    Zcoords=zeros(npoints, npoints);
    
    for r = 1:npoints
        for c = 1:npoints
            Xcoords(r,c) = radius*cos((2*pi*c)/(npoints-1))*sin(pi*r/(npoints-1));
            Ycoords(r,c) = radius*sin((2*pi*c)/(npoints-1))*sin(pi*r/(npoints-1));
            Zcoords(r,c) = radius*cos((pi*r)/(npoints-1));
        end
    end
end

function normal_force = get_force(pos, E_effective, radius)
    
    if pos <= radius
        d = radius - pos;
        k = 1;
    end
    if pos > radius
        d = 0;
        k = 1;
    end
    if pos < 0
        k = d;
        d = radius;
    end

    %in the quasi-static case, the surface integral evaluates to this nice
    %formula. if we induce waves onto the profile, we will need to compute
    %the surface integral manually within this function
    normal_force = (2*k-2/3)*E_effective*sqrt(radius*d^3);
end
