
%This is the folder that we will put the positions into
fileID = fopen('points.txt', 'w');
numballsinstring = 158;%Gives a nice value of the loop thickness d

%We will go around the circle and place particles for to be used as
%strings or rings. First we will create a theta and dtheta which will go around the
%circle and point at each point. Then we will calculate the x and y values
%from the value of theta. We want to make strings which are a tiny bit
%smaller than the cylinder. Cylinder is 8cm in diameter, 4cm in radius, so
%we want the length of the pointer to be 
R = (0.08-0.00425)/2; %With a lil buffer
%After a teensy trig
dth = 2*pi/numballsinstring;
Th = 0:dth:(2*pi-dth);
x = R*cos(Th);
y = R*sin(Th);
A = [x; y];
format long
actuald = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2)
actualThnaught = atan2(y(3) - y(1), x(3) - x(1)) - atan2(y(2) - y(1), x(2) - x(1))
fprintf(fileID, 'create_atoms 1 single %13.9f %13.9f ${zlevel}\n', A);
fclose(fileID);

