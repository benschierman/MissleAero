%% Vortex Latice Method - Missle Aerodynamics Project 1: Due Aug. 20, 2021
close all; clear;

%% Panel, Aspect Ratio, and Length Input
% n=input('How many panels would you like to calculate? \n');
% 
% if isempty(n);
    n=4;
% elseif n==0;
%     print('Please enter a value greater than 0. \n')
%     n=input('How many panels would you like to calculate? \n')
% end 
% and now i do an edit
% ar = input('What is your aspect ratio? \n');
% if isempty(ar);
%     ar=3.55;
%     %ar = 5;
% end
% 
% b1 = input('What is your wingspan? \n');
% if isempty(b1);
    b1 = 1;
% end
% 
% tr = input('What is your taper ratio? \n');
% if isempty(tr);
%     %tr = 1;
    tr = 0.5;
% end
% 
% % cr = input('What is your root chord length');
% % if isempty(cr);
% %     cr=.2;
% % end
% 
% theta = input('What is your sweep angle (in degrees)? \n');
% if isempty (theta);
    theta = 45;
%     theta = -45;
% end
ar = 1;
index = 1;
%% Lets try to optimize
while ar < 100
%% Misc. Calculations

ceq=(b1/ar);
cr = (ceq)*2/(1+tr);
ct = cr*tr;
thetaR=deg2rad(theta);
A=b1/(2*n);


%% Plot coordinates via functions, accounting for taper ratio and sweep

%quarter chord

g = @(x) ((.5*b1*tan(thetaR))/(.5*b1))*x+(1/4)*cr;
% fplot(g, [0 .5*b1], '--r') 
% set (gca, 'YDir', 'reverse')
% hold on 
gg = @(x) -((.5*b1*tan(thetaR))/(.5*b1))*x+(1/4)*cr;
% fplot(gg, [-.5*b1 0], '--r')


% 3-Quarter Chord Line

f = @(x) ((.5*b1*tan(thetaR)+(1/2)*ct-(2/4)*cr)/(.5*b1))*x+(3/4)*cr;
% fplot(f, [0 .5*b1], '--b')
ff = @(x) -((.5*b1*tan(thetaR)+(1/2)*ct-(2/4)*cr)/(.5*b1))*x+(3/4)*cr;
% fplot(ff, [-.5*b1 0], '--b')

% wings

h = @(x) ((g(.5*b1)-(1/4)*ct)/(.5*b1))*x;
hh = @(x) -((g(.5*b1)-(1/4)*ct)/(.5*b1))*x;
% fplot(h, [0 0.5*b1], 'k')
% fplot(hh, [-.5*b1 0], 'k')

j = @(x) ((((g(.5*b1)+(3/4)*ct))-cr)/(.5*b1))*x + cr;
jj= @(x) -((((g(.5*b1)+(3/4)*ct))-cr)/(.5*b1))*x + cr;
% fplot(j, [0 0.5*b1], 'k')
% fplot(jj, [-.5*b1 0], 'k')

rstx = .5*b1;
rsty = g(rstx) - (1/4)*ct;
rsbx = rstx;
rsby = g(rstx) + (3/4)*ct;
% plot([rstx rsbx], [rsty rsby], 'k')
% plot([-rstx -rsbx], [rsty rsby], 'k')
% xlabel('y', 'FontSize', 20);
% ylabel('x', 'FontSize', 20)


%% Finding individual points and plotting
for i=1:n
    y1n(i) = (i-1)*(.5*b1)/n;
    y2n(i) = (i-1)*(.5*b1)/n + (.5*b1)/n;
    ym(i) = (i-1)*(.5*b1)/n + (.5*b1)/(2*n);
end
    x1n = g(y1n);
    x2n = g(y2n);
    xm = f(ym);
    
% plot(ym, xm, 'or')
% plot(y1n, x1n, '*b')
% plot(y2n, x2n, '*c')
% plot(-y1n, x1n, '*b')
% plot(-y2n, x2n, '*c')
% plot(-ym, xm, 'or')
% 
% legend('Quarter Chord Line', '', 'Three-Quarter Chord Line', '', 'Wings', '', '', '', '', '', '', 'XM1, YM1', 'XM2, YM2', ...
%     '', '', 'Control Points', 'FontSize', 20, 'Location', 'North')

%% Calculating ws

for i=1:n
    for ii=1:n 
        a(ii) = 1/((xm(ii)-x1n(i))*(ym(ii)-(y2n(i)))-(xm(ii)-x2n(i))*(ym(ii)-y1n(i)));
        b(ii) = (((x2n(i)-x1n(i))*(xm(ii)-x1n(i))+(y2n(i)-y1n(i))*((ym(ii)-y1n(i)))))/sqrt((xm(ii)-x1n(i))^2+(ym(ii)-y1n(i))^2);
        c(ii) = ((x2n(i)-x1n(i))*(xm(ii)-x2n(i))+(y2n(i)-y1n(i))*((ym(ii)-y2n(i))))/(sqrt((xm(ii)-x2n(i))^2+(ym(ii)-y2n(i))^2));
        d(ii) = (1/(y1n(i)-ym(ii)))*(1+((xm(ii)-x1n(i))/(sqrt((xm(ii)-x1n(i))^2+(ym(ii)-y1n(i))^2))));
        e(ii) = (1/(y2n(i)-ym(ii)))*(1+((xm(ii)-x2n(i))/(sqrt((xm(ii)-x2n(i))^2+(ym(ii)-y2n(i))^2))));
        ws(ii,i) = a(ii)*(b(ii)-c(ii))+d(ii)-e(ii);
    end
end

%% Switching y values for port side
y1n=-y1n;
y2n=-y2n;

%% Calculating wp
for i=1:n
    for ii=1:n 
        a(ii) = 1/((xm(ii)-x1n(i))*(ym(ii)-(y2n(i)))-(xm(ii)-x2n(i))*(ym(ii)-y1n(i)));
        b(ii) = (((x2n(i)-x1n(i))*(xm(ii)-x1n(i))+(y2n(i)-y1n(i))*((ym(ii)-y1n(i)))))/sqrt((xm(ii)-x1n(i))^2+(ym(ii)-y1n(i))^2);
        c(ii) = ((x2n(i)-x1n(i))*(xm(ii)-x2n(i))+(y2n(i)-y1n(i))*((ym(ii)-y2n(i))))/(sqrt((xm(ii)-x2n(i))^2+(ym(ii)-y2n(i))^2));
        d(ii) = (1/(y1n(i)-ym(ii)))*(1+((xm(ii)-x1n(i))/(sqrt((xm(ii)-x1n(i))^2+(ym(ii)-y1n(i))^2))));
        e(ii) = (1/(y2n(i)-ym(ii)))*(1+((xm(ii)-x2n(i))/(sqrt((xm(ii)-x2n(i))^2+(ym(ii)-y2n(i))^2))));
        wp(ii,i) = a(ii)*(b(ii)-c(ii))+d(ii)-e(ii);
    end
end

%% Flip the sign? I guess it makes sense that it's negative since you sum them?
ws;
wp=-wp;
w=ws+wp;

%% Determine the strength of each of the vortices (gamma) by taking the inverse. I for some reason have to make it negative to get the right answer.
gamma=-sum(inv(w), 2); %true gamma is this multiplied by 4*pi*b*U_inf*alpha

%% Calculate lift (l)
l=sum(gamma); %true lift is this multiplied by p_inf*U_inf^2*b^2*pi*alpha

%% Lift coefficient (cl)
l1=l*2*ar;
cl=l1*pi; %this is PER radians
cld=(cl/180)*pi; %this is PER degree

cldi(index) = cld;
ar = ar+0.25
index = index +1;

end
%% Plot lift coefficient
% figure(2);
% liftcoef = @(x) cld*x;
% fplot(liftcoef, [0 10])
% dots = [2:2:10];
% dotsy = liftcoef(dots);
% 
% hold on 
% axis([0 12 0 .6])
% 
% plot(dots, dotsy, 'ok')
% ylabel('C_L', 'FontSize', 24)
% xlabel('\alpha, deg', 'FontSize', 24)
% 
% legend('C_L/\alpha, deg', '', 'FontSize', 24)


