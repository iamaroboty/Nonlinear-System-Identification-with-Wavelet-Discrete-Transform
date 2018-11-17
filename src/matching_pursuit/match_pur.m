%% matching pursuit basis comparison 


% possible non-linearities 

load cuspamax;
lstcpt = {'poly'};
mpdict = wmpdictionary(length(cuspamax),'LstCpt',lstcpt);
[yfit,r,coeff,iopt,qual] = wmpalg('OMP',cuspamax,...
    mpdict,'wmpcfs',0.8, 'typeplot','stepwise','stepplot',10);