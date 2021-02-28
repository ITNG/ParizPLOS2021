% Program for loading txt files and changing them to mat
ensembles=4;
RHO=cell(1,ensembles);
% RR=cell(1,ensembles);
for ens=1:ensembles
    [dee,dI]=g_extractor(ens);
    fname0=['rho_deeout_',num2str(dee),'_ens',num2str(ens),'_dI',num2str(dI),'.txt'];
    if exist(fname0,'file')
        rho=load(fname0);
        rho=logical(rho');
        display('rho loaded')
    end
    RHO{ens}=rho;
end
fname=['RES_de',num2str(dee),'_dIn_',num2str(dI),'.mat'];
save(fname,'RHO')
