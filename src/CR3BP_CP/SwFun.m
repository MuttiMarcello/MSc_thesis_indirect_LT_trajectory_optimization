function [S] = SwFun(tt,zz,prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% throttle switching function for TO (isFO=0) and EO/FO problems (isFO=1)

    isFO=prob.isFO;
    
    S=zeros(length(tt),1);
    
    for i=1:length(tt)
        
        if isscalar(tt)
            m=zz(7);
            llv=zz(11:13);
            lm=zz(14);
        else
            m=zz(i,7);
            llv=zz(i,11:13);
            lm=zz(i,14);
        end

        c=prob.param(2);
    
        lv=norm(llv);
        
        S(i)=-lv.*c./m-lm+isFO;
    
    end

end

