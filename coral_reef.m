%% Simulating growth and drowning of coral reefs - week 3
% sinusoidal sea level, downgoing linear slab, growing coral
% JSB 1/30/16

clear all
figure(3)
clf

%% Initialize

dx = 100; %  x step m
xmin = 40*1000; % x min distance km
xmax = 100 * 1000; % x max distance km
x = xmin:dx:xmax;% distance array

dz = 1; % z step m
zmin = -1500; % depth min m
zmax = 500; % depth max m
z = zmin:dz:zmax; % depth array

dt = 100; % time step year
tmax = 100*1000; % max time ka
t = 0:dt:tmax; % time array
time = 0;

meansl = 0; % mean sea level

A = 60; % sea level amplitude m
P = 25 * 1000; % sea level change period years

% con_rate = .1; % convergence rate 10cm/year
sub_rate = .008; % subsidence rate 8mm/year

Gm = .01; % max upward growth rate (10 - 15 mm/year)
k = .1; % extinction coeficient (.04 - 0.16 1/m)
I0 = (2000*(10^-6))*(365*24*3600); % surface light intensity (2000 - 2250 microJ/m2*year)
Ik = (450*(10^-6))*(365*24*3600); % saturating light intensity (50-450- microJ/m2*year)

imax = length(t);
%% Run

plate = (-0.05*x)+ 4000; % initial plate geometry
coral = plate; % initial coral geometry - grow on plate

for i = 1:imax  % time loop
sealevel = (0*x)+meansl+A*sin(2*pi*(t(i)/P)); % absolute change of sea level    
z = sealevel-coral; % change in z due to sea level
G = Gm*tanh((I0*(exp(-k*z)))/Ik); % growth rate from Galewsky calculated with changing z
G(coral>=sealevel)=0; % can't grow above water
plate = plate-(dt*sub_rate); % subsiding plate over time
coral = coral+(G*dt)-(dt*sub_rate); % coral growth over time
time = time+dt;

%sealevelrate = t(i)*(1/P)*A*2*pi*cos(2*pi*(t/P)); % rate of change of sealevel

figure(3)
plot(x,sealevel,'-.')
hold on
%plot(G,z,'g','linewidth',1)
plot(x,coral,'m','linewidth',2)
hold on
plot(x,plate,'k','linewidth',2)
hold off
   xlabel('Distance (m)','fontname','arial','fontsize',21)
   ylabel('Depth (m)','fontname','arial','fontsize',21)
   set(gca,'fontsize',18,'fontname','arial')
   title(['Coral evolution after ',num2str(time),' years'])
   
axis([xmin xmax zmin zmax])

   pause(.05)
   
end
