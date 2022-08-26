% load saved figures
fig1   = hgload('/path');
fig2   = hgload('/path');
%  Prepare 'subplot'
figure
h1=subplot(1,1,1);
h2=subplot(1,1,1);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
hold on
copyobj(allchild(get(fig2,'CurrentAxes')),h1)