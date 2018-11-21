%% matching pursuit basis comparison 
clear all
close all


% possible non-linearities 

load cuspamax;
x = gdist(0.1,-100:100);
 x = tanh(-100:100);
lstcpt = {'dct'};
mpdict = wmpdictionary(length(x),'LstCpt',lstcpt);
[yfit,r,coeff,iopt,qual] = wmpalg('OMP',x,...
    mpdict,'wmpcfs',0.8, 'typeplot','movie','stepplot',10);


% load cuspamax;
% x = tanh(0:100);
% x = tube(-100:100, 0.7);
% x = gdist(0.1,-100:100);
% mpdict = wmpdictionary(length(x),'LstCpt',...
%    {'dct'});
% yfit = wmpalg('OMP',x,mpdict);
% plot(x,'k'); hold on;
% plot(yfit,'linewidth',2); legend('Original Signal',...
%     'Matching Pursuit');


% [mpdict,~,~,longs] = wmpdictionary(50,'lstcpt',{{'haar',2}});
% 
% for nn = 1:size(mpdict,2)
%     if (nn<=longs{1}(1))
%         plot(mpdict(:,nn),'k','linewidth',2)
%         grid on
%         xlabel('Translation')
%         title('Haar Scaling Function - Level 2')
%     elseif (nn>longs{1}(1) && nn<=longs{1}(1)+longs{1}(2))
%         plot(mpdict(:,nn),'r','linewidth',2)
%         grid on
%         xlabel('Translation')
%         title('Haar Wavelet - Level 2')
%     else
%         plot(mpdict(:,nn),'b','linewidth',2)
%         grid on
%         xlabel('Translation')
%         title('Haar Wavelet - Level 1')
%     end
%     pause(0.2)
% end