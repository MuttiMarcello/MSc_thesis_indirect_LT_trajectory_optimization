function [prob] = TO_sol(prob)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

% single TO solutions search

    fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-14,'FunctionTolerance',1e-14,'MaxIterations',2e2);

    TU=3.75676967e5; % [s]

    it=1;

    t_wo=prob(it).tw(1);

    prob(it).t0=t_wo;

    N=25; % init solution monotonic search discretization grid

    ex_flag=0;
    f=1;
    atmp=1; % number of rundown attempts
    tf_gv=linspace(1,15,N)*86400/TU; % rundown guesses [1-15 days]

    wb1=waitbar(0,'Generating first TO solution');

    lltf_TO=zeros(8,1);

    while ex_flag<=0

        tf_g=tf_gv(f);
        % tf_g=7.8549*86400/TU;

        lltf_g=[ACT(prob(it)); tf_g];
        % lltf_g(1:7)=[ -28.9365
        %               -74.4595
        %                -0.5211
        %                 0.2184
        %                -0.1160
        %                -0.0001
        %                 0.9977];

        [lltf_TO,~,ex_flag]=fsolve(@(llt) TO_ZFP(llt,prob(it)),lltf_g,fsopt);

        

        wb1=waitbar(f/(atmp*N),wb1,'Generating first TO solution');

        f=f+1;

        if f>atmp*N
            atmp=atmp+1;
            tf_gv=2*pi*linspace(1,2,atmp*N);
            f=1;
        end

    end

    [~,~,prob(it)]=TO_ZFP(lltf_TO,prob(it));
    prob(it)=DispRes(prob(it),1);

    close(wb1);

end