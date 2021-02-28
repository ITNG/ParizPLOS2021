function A=create_A(num_net,N,Nep,inter_c,neighbor_c,con)
Ne=Nep*N;
Ni=N-Ne;
A=zeros(num_net*N,'logical');

for nn=0:num_net-1
    for mm=0:num_net-1
        if con(nn+1,mm+1)
            AA=zeros(N);
            for ii=1:N
                if nn==mm
                    p_con=inter_c;
                elseif nn~=mm
                    p_con=neighbor_c;
                end
                x=randi([1 Ne],1,fix(p_con*Ne));
                AA(ii,x)=1;
                if Ni~=0 && nn==mm
                    y=randi([Ne+1 N],1,fix(p_con*Ni));
                    AA(ii,y)=1;
                end
            end
            A(nn*Ne+1:nn*Ne+Ne,mm*Ne+1:mm*Ne+Ne)=AA(1:Ne,1:Ne);
            if Ni~=0 && nn==mm
                A(num_net*Ne+nn*Ni+1:num_net*Ne+nn*Ni+Ni,nn*Ne+1:nn*Ne+Ne)=AA(Ne+1:N,1:Ne);
                A(nn*Ne+1:nn*Ne+Ne,num_net*Ne+mm*Ni+1:num_net*Ne+mm*Ni+Ni)=AA(1:Ne,Ne+1:N);
                A(num_net*Ne+mm*Ni+1:num_net*Ne+mm*Ni+Ni,num_net*Ne+nn*Ni+1:num_net*Ne+nn*Ni+Ni)=AA(Ne+1:N,Ne+1:N);
            end
        end
    end
end
A=A.*~eye(num_net*N);
