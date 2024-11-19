% MATLAB file 26/05/23
% Adapted from Figure 1 - Izhikevich E.M. (2004) 
% Which Model to Use For Cortical Spiking Neurons? 
%% Limpieza general
clc;clear;close

%% (RS) regular spiking %%%%%%%%%%%%%%%%%%%%%%
a=0.02; b=0.2; c=-65; d=8;
V=-70; u=b*V;
VV=[]; uu=[];
tau = 0.25; tspan = 0:tau:100;
T1=tspan(end)/10;
for t=tspan
if (t>T1) 
I=10;
else
I=0;
end
V = V + tau*(0.04*V^2+5*V+140-u+I);
u = u + tau*a*(b*V-u);
if V > 30
VV(end+1)=30;
V = c;
u = u + d;
else
VV(end+1)=V;
end
uu(end+1)=u;
end
figure();
plot(tspan,VV,[0 T1 T1 max(tspan)],-90+[0 0 10 10],'Color',[255 80 80]/255, 'LineWidth',2);
hold on
plot([0 T1 T1 max(tspan)],-90+[0 0 10 10],'Color',[0 0 0],'LineWidth',2);
hold on
% Scale bar
plot([max(tspan)-20 max(tspan)], [-90 -90], 'k', 'LineWidth', 3) % 20 ms
text(max(tspan)-20, -95, ' 20 ms','FontSize', 10, 'VerticalAlignment', 'middle');
text(tspan(1), VV(1)+5, 'V(t)', 'FontSize', 12, 'Color', [255 80 80]/255);
text(0, -85, 'I(t)', 'FontSize', 12, 'Color', [0 0 0]);
hold off
axis off;
% title('Regular spiking (RS)','FontSize',28);
saveas(gcf, 'NeuronRS.svg');

%% (CH) chattering %%%%%%%%%%%%%%%%%%%%%%
a=0.02; b=0.2; c=-50; d=2;
V=-70; u=b*V;
VV=[]; uu=[];
tau = 0.25; tspan = 0:tau:100;
T1=tspan(end)/10;
for t=tspan
if (t>T1) 
I=10;
else
I=0;
end
V = V + tau*(0.04*V^2+5*V+140-u+I);
u = u + tau*a*(b*V-u);
if V > 30
VV(end+1)=30;
V = c;
u = u + d;
else
VV(end+1)=V;
end
uu(end+1)=u;
end
figure();
plot(tspan,VV,[0 T1 T1 max(tspan)],-90+[0 0 10 10],'Color',[192 0 0]/255,'LineWidth',2);
hold on
plot([0 T1 T1 max(tspan)],-90+[0 0 10 10],'Color',[0 0 0],'LineWidth',2);
hold on
% Scale bar
plot([max(tspan)-20 max(tspan)], [-90 -90], 'k', 'LineWidth', 3) % 20 ms
text(max(tspan)-20, -95, ' 20 ms','FontSize', 10, 'VerticalAlignment', 'middle');
text(tspan(1), VV(1)+5, 'V(t)', 'FontSize', 12, 'Color', [192 0 0]/255);
text(0, -85, 'I(t)', 'FontSize', 12, 'Color', [0 0 0]);
hold off
axis off;
% title('Chattering (CH)','FontSize',28);
saveas(gcf, 'NeuronCH.svg');

%% (FS) fast spiking %%%%%%%%%%%%%%%%%%%%%%
a=0.1; b=0.2; c=-65; d=2;
V=-70; u=b*V;
VV=[]; uu=[];
tau = 0.25; tspan = 0:tau:100;
T1=tspan(end)/10;
for t=tspan
if (t>T1) 
I=10;
else
I=0;
end
V = V + tau*(0.04*V^2+5*V+140-u+I);
u = u + tau*a*(b*V-u);
if V > 30
VV(end+1)=30;
V = c;
u = u + d;
else
VV(end+1)=V;
end
uu(end+1)=u;
end
figure();
plot(tspan,VV,[0 T1 T1 max(tspan)],-90+[0 0 10 10],'Color',[102 153 255]/255,'LineWidth',2);
hold on
plot([0 T1 T1 max(tspan)],-90+[0 0 10 10],'Color',[0 0 0],'LineWidth',2);
hold on
% Scale bar
plot([max(tspan)-20 max(tspan)], [-90 -90], 'k', 'LineWidth', 3) % 20 ms
text(max(tspan)-20, -95, ' 20 ms','FontSize', 10, 'VerticalAlignment', 'middle');
text(max(tspan)-20, -95, ' 20 ms','FontSize', 10, 'VerticalAlignment', 'middle');
text(tspan(1), VV(1)+5, 'V(t)', 'FontSize', 12, 'Color', [102 153 255]/255);
text(0, -85, 'I(t)', 'FontSize', 12, 'Color', [0 0 0]);
hold off
axis off;
% title('Fast spiking (FS)','FontSize',28);
saveas(gcf, 'NeuronFS.svg');

%% (AC) accomodation %%%%%%%%%%%%%%%%%%%%%%
a=0.02; b=1; c=-55; d=4;
V=-65; u=-16;
VV=[]; uu=[]; II=[];
tau = 0.5; tspan = 0:tau:400;
for t=tspan
if (t < 200)
I=t/25;
elseif t < 300
I=0;
elseif t < 312.5
I=(t-300)/12.5*4;
else
I=0;
end
V = V + tau*(0.04*V^2+5*V+140-u+I);
u = u + tau*a*(b*(V+65));
if V > 30
VV(end+1)=30;
V = c;
u = u + d;
else
VV(end+1)=V;
end
uu(end+1)=u;
II(end+1)=I;
end
figure()
plot(tspan,VV,'Color',[51 51 255]/255,'LineWidth',2);
hold on
plot(tspan,II*1.5-90,'Color',[0 0 0],'LineWidth',2);
hold on
% Scale bar
plot([max(tspan)-20 max(tspan)], [-90 -90], 'k', 'LineWidth', 3) % 20 ms
text(max(tspan)-20, -95, ' 20 ms','FontSize', 10, 'VerticalAlignment', 'middle');
text(tspan(1), VV(1)+5, 'V(t)', 'FontSize', 12, 'Color', [51 51 255]/255);
text(0, -85, 'I(t)', 'FontSize', 12, 'Color', [0 0 0]);
hold off
axis([0 max(tspan) -90 30])
axis off;
% title('Accomodating (AC)','FontSize',28);
saveas(gcf, 'NeuronAC.svg');
