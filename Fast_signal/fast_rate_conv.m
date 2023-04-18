dh=0.01;Ne=80; 
start_time=1;
end_time=length(rho);
range=length(rho);
% sigma=50;
sigma_w=sigma/dh;
r_sum=Ne.\[double(sum(rho(1:Ne,start_time:end_time),1));double(sum(rho((numnet-1)*Ne+1:numnet*Ne,start_time:end_time),1))];
dhsigma2=2*sigma_w*sigma_w;
B=1000/(sqrt(2*pi)*sigma);
% tic
interval=10*sigma_w;
exp_factor=B*exp(-((-interval/2:interval/2).^2)/dhsigma2);
rr1=conv(r_sum(1,:),exp_factor,'valid');
rr2=conv(r_sum(2,:),exp_factor,'valid');
rr=[rr1;rr2];
% toc
%% 
% figure;
% yyaxis left
% plot(rr(1,:),'r');hold on;plot(rr(2,:),'b');
% x1=(end_time-size(rr,2))/2;
% yyaxis right
% plot(Isig(x1+1:end-x1)+mean(rr(1,5e4:40e4)),'k')
% % findpeaks(rr(1,:),'MinPeakProminence',1,'Annotate','extents');
% % findpeaks(Isig(1:end)+mean(rr(1,5e4:10e4)),'MinPeakProminence',1,'Annotate','extents');
% grid minor
% % ylim([70 75])
% xlim([0 450e3])