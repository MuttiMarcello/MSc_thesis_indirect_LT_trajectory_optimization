function [S] = SwFun(tt,zz,isFO)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP constant parameters

% throttle switching function for TO (isFO=0) and EO/FO problems (isFO=1)
    
    S=zeros(length(tt),1);
    
    for i=1:length(tt)
        
        if isscalar(tt)
            r=norm(zz(1:3));
            m=zz(7);
            llv=zz(11:13);
            lm=zz(14);
        else
            r=norm(zz(i,1:3));
            m=zz(i,7);
            llv=zz(i,11:13);
            lm=zz(i,14);
        end

        Tc=MARGO_param(r);
    
        c=Tc(2);
    
        lv=norm(llv);
        
        S(i)=-lv.*c./m-lm+isFO;
    
    end

end

