npoints = 50;

dx = 1/npoints;
radius = 1;
tfinal = 100;
nsteps = 10000;
dt = tfinal/nsteps;
g = 10;
tau = 1;
m = .1;
U = 0; %elastic potential
pos = 5;
vel = -5;
field = 5;

E_ball = 10^5; %Young's modulus of the ball's material (rubber)
E_wall = 23; %Young's modulus of the wall's material (concrete)
G_ball = .0006; %Shear modulus of the ball's material (rubber)
G_wall = .8; %Shear modulus of the wall's material

E_effective = 1/((1-G_ball^2)/E_ball + (1-G_wall^2)/E_wall);

delta_k = E_effective*dx; %stiffness of spring over area element

KE = 0;

[x,y,z] = create_ball(ceil(npoints/2), radius);
z = z + pos;

x_hs = linspace(-field, field, npoints);
y_hs = linspace(-field, field, npoints);
z_hs = zeros(npoints);

H_hs_prev = zeros(npoints);

for t = 1:nsteps
    clf
    
    KE = .5*m*vel^2;
    normal_force = get_force(pos, E_effective, radius);
    a = normal_force/m - g;
    vel = vel + a*dt;
    %equations needed for inelastic collisions (in which KE is nonconstant)
    %U = U + normal_force*dt;
    %KE = KE - normal_force*dt;
    %vel = sign*sqrt(2*abs(KE)/m) - g*dt;
    pos = pos + vel*dt;
    z = z + vel*dt;
    
    H_hs = profile_ehs(x_hs, y_hs, pos, radius, npoints);
    %v term in case we want to model extra waves...
    %v_hs = (H_hs-H_hs_prev);
    
    surfl(x,y,z); hold all
    surfl(x_hs, y_hs, H_hs); hold all
        axis([-field,field,-field,field,0,2*field])
    %if mod(t, 2) == 0
    drawnow
    %end
end    

%create the sphere
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

%compute the normal force
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

%get the profile of the elastic half space
function ehs = profile_ehs(x_hs, y_hs, pos, radius, npoints)
    ehs = zeros(npoints);
    
    if pos > radius
        return
    else
        a = sqrt(radius^2-(pos)^2);
        r = zeros(npoints);
        
        for i = 1:npoints
            for j = 1:npoints   
                r(i,j)= sqrt(x_hs(i)^2+y_hs(j)^2);
                if r(i,j) < a
                    ehs(i,j) = sqrt(radius^2-(pos)^2);
                elseif r(i,j) < radius
                    ehs(i,j) = (2/pi)*((atanh(a/radius)*(a*asin(a/r(i,j))+sqrt(r(i,j)^2-a^2)))...
                    -radius*asin(a/r(i,j)) + sqrt(radius^2-r(i,j)^2)...
                    *atan((a*sqrt(radius^2-r(i,j)^2))/(radius*sqrt(r(i,j)^2-a^2))));
                elseif r(i,j) >= radius
                    ehs(i,j) = (2/pi)*((atanh(a/radius)*(a*asin(a/r(i,j))+sqrt(r(i,j)^2-a^2)))...
                    -radius*asin(a/r(i,j)));
                end
            end
        end
    end
end