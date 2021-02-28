function J=create_synapse(A,num_net,N,Nep,EEJ,EIJ,IEJ,IIJ,outEEJ,con)
% global A N
Ne=Nep*N;
% Ni=Nip*N;
J=single(zeros(num_net*N));

for nn=0:num_net-1
    for mm=0:num_net-1
        if con(nn+1,mm+1)
            for ii=1:N
                if nn==mm
                    d=EEJ;
                elseif nn~=mm
                    d=EIJ;
                end
                J(nn*Ne+1:nn*Ne+Ne,mm*Ne+1:mm*Ne+Ne)=d;
            end
        end
    end
end
J(:,num_net*Ne+1:end)=EEJ;
J(num_net*Ne+1:end,1:num_net*Ne)=EEJ;
J=J.*A;
