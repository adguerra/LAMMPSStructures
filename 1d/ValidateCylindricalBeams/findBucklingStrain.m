
%I might suggest putting the input file into a folder before you run it, to contain the output. Youll name that folder here
filepath = './';

ng = 10:2:80; %The the set of lengths in number of particles that you simulated 
d = 0.001; %Diameter of particles
l = ng*d; %Length of the beams
bs = zeros(1, length(l)); %Buckling strain
ar = (2*l./d); %aspect ratio
for i = 1:length(l) 
    filename = strcat(filepath, num2str(l(i)), '_time_dl_xofcenter_BE_SE.out'); %Pick the file
    fid= fopen(filename); %open it 
    A = textscan(fid, '%f %f %f %f %f', 'headerlines', 1); %Grab the info out of it 
    fclose(fid); %close it 
    xcenter = A{3}; %take the x position of the center grain
    strain = A{2}./(l(i)); %Calculate the strain
    for j = 1:length(xcenter)
        if abs(xcenter(j)) > 10^-6 %If the center grain moves sideways ...
            continue %...the beam has buckled ...
        end   
        remember = j; %...so remember the frame right before it buckled...
    end
    bs(i) = strain(remember); %... that is your buckling strain
end
%I will also plot the Euler buckling strain, which will miss both our data and Comsols by a lot
x = 20:200;
eps = ((x.^2)/(2*pi^2)).^(-1);
%I manually input these values from a COMSOL simulation
arcomsol = [20 30 40 50 60 70 80 90 100 120 140];
Bstresscomsol = [22048 10244 5853 3769 2626 1934 1482 1172 949.4 659.1 485];
Bscomsol = Bstresscomsol./(960000);
%plot everything
figure();
hold on;
plot(ar, bs, 'o', 'MarkerSize',10);
plot(x, eps);
plot(arcomsol, Bscomsol, 'o', 'MarkerSize',10);
xlim([0 160]);
ylim([0 .025]);
title('Comparrison between Comsol (yellow), Euler (orange), and our beams (blue)')
xlabel('l/r') 
ylabel('Critical Buckling Strain') 

figure();
loglog(ar, bs, 'o', 'MarkerSize',10);
hold on;
loglog(x, eps);
loglog(arcomsol, Bscomsol, 'o', 'MarkerSize',10);
title('Comparrison between Comsol (yellow), Euler (orange), and our beams (blue)')
xlabel('l/r') 
ylabel('Critical Buckling Strain') 
%csvwrite('lammps_aspect_critBStrain.csv', [ar' bs']);
%csvwrite('comsol_aspect_critBStrain.csv', [arcomsol' Bscomsol']);

