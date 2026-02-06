function [H] = Hamil(tt,zz,u,prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP variable parameters

% Hamiltonian from adimensional integrated state

    H=zeros(length(tt),1);

    for i=1:length(H)

        r=norm(zz(i,1:3));
        [~,~,Sp]=MARGO_param(r);
    
        if Sp<prob.Plim(2)
            Ptype='med';
        elseif Sp>=prob.Plim(2)
            Ptype='max';
        end

        ep=prob.epsilon;

        if u(i)==1
            utype='on';
        elseif u(i)==0
            utype='off';
        else
            utype='med';
        end

        if prob.isFO==1

            ff=FO_2BP_SEP(tt(i),zz(i,:).',Ptype,utype,ep);

        else

            ff=TO_2BP_SEP(tt(i),zz(i,:).',Ptype);

        end

        Tc=MARGO_param(r);
        T=Tc(1);
        c=Tc(2);

        ffx=ff(1:7);
        ll=zz(i,8:14).';
        H(i)=dot(ll,ffx)+prob.isFO*(T/c*(u(i)-u(i)*ep*(1-u(i))));
    end
    
end

