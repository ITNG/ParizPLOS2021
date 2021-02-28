function [dee, dI]=g_extractor(x)
xname=ls(['*ens',num2str(x),'*.txt']);
xname(isspace(xname))=[];
dash_place=find(xname=='_');
xdot=find(xname=='.');
dee=str2double(xname(dash_place(2)+1:dash_place(3)-1));
dI=str2double(xname(dash_place(4)+3:xdot-1));
