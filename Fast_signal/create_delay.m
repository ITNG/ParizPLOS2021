function delay=create_delay(A,num_net,N,Nep,in_delay,out_delay,con)
% global A N
Ne=Nep*N;
% Ni=Nip*N;
delay=single(zeros(num_net*N));

for nn=0:num_net-1
    for mm=0:num_net-1
        if con(nn+1,mm+1)
            for ii=1:N
                if nn==mm
                    d=in_delay;
                elseif nn~=mm
                    d=out_delay;
                end
                delay(nn*Ne+1:nn*Ne+Ne,mm*Ne+1:mm*Ne+Ne)=d;
            end
        end
    end
end
delay(:,num_net*Ne+1:end)=in_delay;
delay(num_net*Ne+1:end,1:num_net*Ne)=in_delay;
delay=delay.*A;
