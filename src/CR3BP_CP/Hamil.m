function [H] = Hamil(tt,zz,u,prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% Hamiltonian from adimensional integrated state

    H=zeros(length(tt),1);

    for i=1:length(H)

        ep=prob.epsilon;


        if prob.isFO==1

            ff=FO_CR3BP_cp(tt(i),zz(i,:).',prob);

        else

            ff=TO_CR3BP_cp(tt(i),zz(i,:).',prob);

        end

        T=prob.param(1);
        c=prob.param(2);

        ffx=ff(1:7);
        ll=zz(i,8:14).';
        H(i)=dot(ll,ffx)+(1-prob.isFO)+prob.isFO*(T/c*(u(i)-u(i)*ep*(1-u(i))));
    end
    
end

