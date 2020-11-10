%% AA228 Final Project Martin Kamme
close all
clear all
clc

muearth = 398600;
rearth = 6378;

%% Setup Chaser (A) and Target (B)
incA = 0;
eccA = 0;
RAANA = deg2rad(115.4157);        % radians
perA = deg2rad(75.9974);          % radians
thetaA  = deg2rad(190.4478);      % radians
TA = 23.94*60*60;                 % period in seconds
nA = ((2*pi)/TA);                 % mean motion of the target
aA = (muearth/(nA^2))^(1/3);      % km
hA = sqrt(aA*muearth*(1-eccA^2));
[rinitialA,vinitialA] = coes2rv(hA,incA,RAANA,eccA,perA,thetaA,muearth);

incB = 0;
eccB = 0;
RAANB = deg2rad(115.4157);        % radians
perB = deg2rad(75.9974);          % radians
TB = 23.94*60*60;                 % period in seconds
nB = ((2*pi)/TB);                 % mean motion of the target
aB = (muearth/(nB^2))^(1/3);      % km
hB = sqrt(aB*muearth*(1-eccB^2));

rHill = [0;10;0];                 % km
vHill = [0;0;0];                  % km/s


%% Parameters of planning
% Tcoast;
% numburns;
% dvBurn;
% targetSeparation;

%% Start with making some observations
% make the observations by generating output with the CWHpropagator and
% adding gaussian noise to the output
clc

% Burn calc
dVY = (rHill(2) - 5)*(nA/(6*pi));
vHill = vHill + [0;dVY;0];

t = linspace(1,TA,100000);

% Standard Deviations of noisy data
stdX = .1; % km
stdY = .1; % km

% Generate the true path and noisy observations
for i = 1:length(t)
    [rHillplot(:,i),vHillplot(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
    rHillplotNoisy(1,i) = rHillplot(1,i) + stdX*randn(size(rHillplot(1,i)));
    rHillplotNoisy(2,i) = rHillplot(2,i) + stdY*randn(size(rHillplot(2,i)));
    %rHillplotNoisy(3,i) = rHillplot(3,i) + stdZ*randn(size(rHillplot(3,i)));
    %vHillplotNoisy(1,i) = vHillplot(1,i) + stdX*randn(size(vHillplot(1,i)));
    %vHillplotNoisy(2,i) = vHillplot(2,i) + stdY*randn(size(vHillplot(2,i)));
    %vHillplotNoisy(3,i) = vHillplot(3,i) + stdZ*randn(size(vHillplot(3,i)));
    
    if i ~= 1
        dataOut(i,:) = [t(i)-t(i-1) , rHillplotNoisy(:,i)']; %, vHillplotNoisy(:,i)'];
    else
        dataOut(i,:) = [0 , rHillplotNoisy(:,i)'];% , vHillplotNoisy(:,i)'];
    end
    
end

% Write the data to a csv
filename = 'HopObservations.xlsx';
recycle on
delete(filename)
writematrix(dataOut,filename)


%% Run the LKF
clc
X0 = [rHillplotNoisy(:,1);vHillplot(:,1)];
% uncertainties of initial conditions
xunc = .05; % km
yunc = .05; % km
xdotunc = .001; % km/s
ydotunc = .001; % km/s
P0 = [xunc^2 0 0 0; 0 yunc^2 0 0;0 0 ydotunc^2 0 ;0 0 0 ydotunc^2];

observations = readmatrix('HopObservations.xlsx');

Xout = LKFrelative(X0,P0,observations,nA,stdX);

for i = 1:length(Xout)
    errX(i) = Xout(i,1) - rHillplot(1,i);
    errY(i) = Xout(i,2) - rHillplot(2,i);
    rangeError(i) = sqrt(errX(i)^2 + errY(i)^2);
end

%% Plots
clc
close all

originx = 0;
originy = 0;
line0 = [0 TA/60];
line0y = [0 0];
figure
hold on
plot(rHillplotNoisy(2,:),rHillplotNoisy(1,:),'rx')
plot(Xout(:,2),Xout(:,1),'mo')
plot(rHillplot(2,:),rHillplot(1,:),'b-','LineWidth',2)
plot(originx,originy,'rs','MarkerFaceColor','r','MarkerSize',10,'MarkerEdgeColor','k')
plot(10,0,'bs','MarkerFaceColor','b','MarkerSize',10,'MarkerEdgeColor','k')
plot(5,0,'bs','MarkerFaceColor','g','MarkerSize',10,'MarkerEdgeColor','k')
xlabel('y (km)') % Vbar
ylabel('x (km)') % Rbar
legend('Noisy Data','Kalman Filtered Position','True Path','Target Spacecraft','Chaser Spacecraft Initial Postion','Chaser Spacecraft Final Position','Location','northeast')
xlim([0 15])
ylim([0 3])

figure
hold on
plot(t./60,errX*1000,'mx-')
plot(line0,line0y,'b','LineWidth',2)
xlabel('time (min)')
ylabel('Error (m)')
legend('Error in x','No error')

figure
hold on
plot(t./60,errY*1000,'kx-')
plot(line0,line0y,'b','LineWidth',2)
xlabel('time (min)')
ylabel('Range Error (m)')
legend('Error in y','no error')

figure
hold on 
plot(t./60,rangeError*1000,'mx-')
xlabel('time (min)')
ylabel('Range Error (m)')
legend('Range Error')


%% Functions

% Linear Kalman Filter for Relative Motion
function Xout = LKFrelative(X0,P0,observations,omega,stdX)
% INPUT:
% X0 contains:
% X0(1) = xhill
% X0(2) = yhill
% X0(3) = zhill
% X0(4) = xdothill
% X0(5) = ydothill
% X0(6) = zdothill
%
% P0 contains: a 4x4 covariance matrix

% partial derivative matrix of the observations
H = [1 0 0 0;
    0 1 0 0];         % only making range observations
sigma = stdX;         % standard deviation of the measurement (km)

% noise statistics, what are the uncertainties of the dynamics model?
xnoise = .001;      % km
ynoise = .001;      % km
xdotnoise = 1e-7;   % km/s
ydotnoise = 1e-7;   % km/s
Q = [xnoise^2 0 0 0 ; 0 ynoise^2 0 0 ; 0 0 xdotnoise^2 0 ; 0 0 0 ydotnoise^2];

% measurement noise matrix
R = [sigma^2 0;
    0 sigma^2];

for i = 1:length(observations)
    
    % Transition matrix
    t = observations(i,1);
    PHI =  [4-3*cos(omega*t) , 0 , 1/omega*sin(omega*t) , -2/omega*cos(omega*t)+2/omega;
        6*sin(omega*t)-6*omega*t , 1 , 2/omega*cos(omega*t)-2/omega , 4/omega*sin(omega*t)-3*t ;
        3*omega*sin(omega*t) , 0  , cos(omega*t) , 2*sin(omega*t);
        6*omega*cos(omega*t)-6*omega , 0 , -2*sin(omega*t) , 4*cos(omega*t)-3];
    
    % Predict the state
    Xpredict = PHI*X0;
    
    % Predict the covarianve matrix
    Ppredict = PHI*P0*PHI' + Q;
    
    % Calculate the Kalman Gain
    K = Ppredict*H'*pinv(H*Ppredict*H'+ R);
    
    % Calculate the state
    X0 = Xpredict + K*([observations(i,2:3) ]' - H*Xpredict);
    
    % Calculate the covariance matrix
    P0 = (eye(4)-K*H)*Ppredict;
    
    Xout(i,:) = X0';
end

end


% Coes to R and V
function [R,V] = coes2rv(h,inc,RAAN,e,per,theta,muearth)
% h [km^2/s] Specific angular momentum
% i [rad] Inclination
% RAAN [rad] Right ascension (RA) of the ascending node
% e Eccentricity
% per [rad] Argument of perigee
% theta [rad] True anomaly
% muearth = 398600; Earth’s gravitational parameter [km^3/s^2]

% State Vectors in Perifocal coordinates
rx = h^2/muearth*(1/(1 + e*cos(theta)))*[cos(theta);sin(theta);0];
vx = muearth/h*[-sin(theta); (e +cos(theta));0];

% Direction cosine matrix
DCM = [cos(per), sin(per),0;-sin(per),cos(per),0;0,0,1]*...
    [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)]*...
    [cos(RAAN), sin(RAAN),0;-sin(RAAN),cos(RAAN),0;0,0,1];

% Transformation Matrix
Dcm = inv(DCM);

% ECI R
R = Dcm*rx;

% ECI V
V = Dcm*vx;

end

function [rHill,vHill] = CWHPropagator(rHillInit,vHillInit,omega,t)
% Purpose:
% Take initial position and velocity coordinates in the Hill reference frame
% and propagate them using the Clohessy-Wiltshire Hill Linearize equation
% of motion.
%
% Inputs:
%rHillInit                  [3 x 1]                 Hill Position vector
%                                                   (km) / (m)
%
%vHillInit                  [3 x 1]                 Hill Velocity vector of
%                                                   (km/s) / (m/s)
%
%omega                       double                 Orbital Angular Rate
%                                                   of the target
%                                                   (rad/s)
%                                                   Should be close to
%                                                   circular for linear propagation
%                                                   error to be low.
%
%t                          [1 x N]                 Propagation Time in
%                                                   seconds
%
%
%
%
% Outputs:
%rHill                       [3 x N]                Propagated Hill
%                                                   Position vector (km) /
%                                                   (m/s)
%
%vHill                       [3 x N]                Propagated Hill
%                                                   Velocity vector (km/s)
%                                                   / (m/s)
%
%
% References:
%
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
x0 = rHillInit(1,:); y0 = rHillInit(2,:); z0 = rHillInit(3,:);
x0dot = vHillInit(1,:); y0dot = vHillInit(2,:); z0dot = vHillInit(3,:);
rHill = [(x0dot./omega).*sin(omega.*t)-(3.*x0+2.*y0dot./omega).*cos(omega.*t)+(4.*x0+2.*y0dot./omega)
    (6.*x0+4.*y0dot./omega).*sin(omega.*t)+2.*(x0dot./omega).*cos(omega.*t)-(6.*omega.*x0+3.*y0dot).*t+(y0-2.*x0dot./omega)];
%z0.*cos(omega.*t)+(z0dot./omega).*sin(omega.*t)];
vHill = [x0dot.*cos(omega.*t)+(3.*omega.*x0+2.*y0dot).*sin(omega.*t)
    (6.*omega.*x0 + 4.*y0dot).*cos(omega.*t) - 2.*x0dot.*sin(omega.*t)-(6.*omega.*x0 + 3.*y0dot)];
%-z0.*omega.*sin(omega.*t)+z0dot.*cos(omega.*t)];
end
