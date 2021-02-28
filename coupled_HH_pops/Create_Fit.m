%% Fit: 'untitled fit 1'.

x=DE;
y=DI;
% CC=mean(zlc_rr2_Isig_ub,3);
% CC(:,end)=[];
% CC(end,:)=[];
[xData, yData, zData] = prepareSurfaceData( x, y, CC );

% Set up fittype and options.
ft = 'biharmonicinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'off' );


x1=linspace(min(x),max(x),500)';
y1=linspace(min(y),max(y),500)';
yy=zeros(numel(x1)*numel(y1),1);xx=zeros(numel(x1)*numel(y1),1);
for ii=1:numel(x1)
    xx((ii-1)*numel(x1)+1:ii*numel(x1))=x1(ii)*ones(numel(y1),1);
    yy((ii-1)*numel(y1)+1:ii*numel(y1))=y1;
end
xZ=fitresult([xx,yy]);
Z=zeros(numel(x1));
for ii=1:numel(x1)
    Z(:,ii)=xZ((ii-1)*numel(y1)+1:ii*numel(y1));
end
figure;
imagesc(x1,y1,Z)