npoints = 25;
radius = 1;
tfinal = 100;
nsteps = 10000;
dt = tfinal/nsteps;
g = 10;
tau = 1;
m = 1;

u = zeros(npoints);
v = zeros(npoints);
w = zeros(npoints);


[x,y,z] = create_ball(npoints, radius);
z = z + 10;

for t = 1:nsteps
    clf
    
    w = w - g*dt;
    z = z + w*dt;
    for i = 1:npoints
        for j = 1:npoints   
            if z(i,j) < 0
                z(i,j) = 0;
            end
        end
    end
    
    surfl(x,y,z); hold all
        axis([-5,5,-5,5,0,10])
    drawnow
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