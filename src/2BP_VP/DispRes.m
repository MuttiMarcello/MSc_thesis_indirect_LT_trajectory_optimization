function [prob] = DispRes(prob,output)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% 2BP variable parameters

% Graphical display of trajectory and quantities of interest

    if nargin<2
        output=1;
    end

    % Set the default font size for axes labels, tick labels and legend
    set(0,'DefaultAxesFontSize',15);
    set(0,'DefaultLegendFontSize',15);
    
    % Set the default linewidth
    set(0,'DefaultLineLineWidth',2);
    
    % Set the default text and legend interpreter
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    g0=9.80665e-3;      % [km/s^2]

    m0=prob.m0;
    t0=prob.t0;
    targ=prob.targ;

    Pmin=prob.Plim(1);
    Pmax=prob.Plim(2);

    % dimensional units
    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1
    MU=prob.m0;

    tf_ad=prob.tf_ad;
    tf=prob.tf;

    y0=prob.y0;

    Phi0=eye(length(y0));
    vPhi0=reshape(Phi0,[length(y0)^2,1]);

    if prob.isFO==1

        [tt,zz]=FO_ode87(prob,[0 tf_ad],[y0; vPhi0]);

    else

        [tt,zz]=TO_ode87(prob,[0 tf_ad],[y0; vPhi0]);

    end

    ToF=(tf-t0)/86400;          % time-of-flight [d]
    mf=zz(end,7)*m0;            % final mass [kg]
    mp=(zz(1,7)-zz(end,7))*m0;  % propellant mass [kg]

    ttd=tt*TU/86400; % [d]

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');
    xtf_ad=ADIM([xtf; 1],m0);
    rvtf=xtf_ad(1:6);

    % terminal target error
    df=zz(end,1:6).'-rvtf;

    r=sqrt(sum(zz(:,1:3).^2,2));

    % propulsive parameters profiles
    [TTcc,~,P]=MARGO_param(r);

    TT=(MU*LU)/(TU^2)*TTcc(1,:)*1e6; % [mN]
    II=LU/TU*TTcc(2,:)/g0;           % [s]

    % throttle switching function profiles
    ep=prob.epsilon;
    Se=SwFun(tt,zz,prob.isFO);

    % throttle profile
    if ep>0

        u=1.*(Se<-ep)+(ep-Se)./(2*ep).*(abs(Se)<=ep)+0;

    else
        
        u=zeros(size(Se));
        ut=1.*(Se(1)<ep)+0;
        u(1)=ut;

        for j=2:length(Se)
            if ttd(j)==ttd(j-1) && isapprox(Se(j),ep,AbsoluteTolerance=1e-10)
                u(j)=1-ut;
                ut=u(j);
            else
                u(j)=ut;
            end
        end

    end

    % hamiltonian
    H=Hamil(tt,zz,u,prob);

    prob.mf=mf;
    prob.mp=mp;
    prob.S=Se;
    prob.H=H;
    prob.tt_ad=tt;
    prob.tt=ttd;
    prob.zz=zz;

    if output==1

        plot3D(t0,ttd,zz,targ,u);
    
        fprintf('Departure date: %s (%.1f MJD2000)\n',cspice_et2utc(t0,'C',0),et2MJD2000(t0))
        fprintf('Arrival date: %s (%.1f MJD2000)\n',cspice_et2utc(tf,'C',0),et2MJD2000(tf))
        fprintf('ToF: %.4f days, Final Mass: %.4f kg, Propellant Mass: %.4f kg\n\n',ToF,mf,mp)
    
        fprintf('Position error %.4e km\nVelocity error %.4e km/s\n',norm(df(1:3))*LU, norm(df(4:6))*LU/TU)

        fprintf('Max relative Hamiltonian variation %.4e%%\n\n',abs((max(H)-min(H))/max(H))*100)
    

        %-throttle---------------------------------------------------------
        subplot(3,3,2)
        plot(ttd,u)
        axis tight
        ylim([0 1.1])
        ylabel('$u$')
        grid on
        grid minor

        %-switching-function-----------------------------------------------
        subplot(3,3,5)
        plot(ttd,Se)
        axis tight
        ylabel('$S$')
        grid on
        grid minor
            
        %-mass-------------------------------------------------------------
        subplot(3,3,8)
        plot(ttd,zz(:,7)*m0)
        axis tight
        ylabel('$m\,[kg]$')
        xlabel('$t\,[days]$')
        grid on
        grid minor

        %-input-power------------------------------------------------------
        subplot(3,3,3)
        plot(ttd,P)
        if max(P)>Pmax
            hold on
            plot([ttd(1) ttd(end)],[Pmax,Pmax],'--r','LineWidth',1)
        end
        if min(P)<Pmin
            hold on
            plot([ttd(1) ttd(end)],[Pmin,Pmin],'--r','LineWidth',1)
        end
        axis tight
        ylabel('$S_{p}\,[W]$')
        grid on
        grid minor

        %-thrust-----------------------------------------------------------
        subplot(3,3,6)
        plot(ttd,II)
        axis tight
        ylabel('$I_{sp}\,[s]$')
        grid on
        grid minor
    
        %-specific-impulse-------------------------------------------------
        subplot(3,3,9)
        plot(ttd,TT)
        axis tight
        ylabel('$T_{max}\,[mN]$')
        xlabel('$t\,[days]$')
        grid on
        grid minor
        
    end



    % User plot setting removal
    set(0,'DefaultAxesFontSize','remove');
    set(0,'DefaultLegendFontSize','remove');
    set(0,'DefaultLineLineWidth','remove');
    set(groot,'defaultTextInterpreter','remove');
    set(groot,'defaultLegendInterpreter','remove');
end

function plot3D(t0,tt,zz,targ,u)
    
    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    ttd=tt.'.*86400+t0;

    xx_SEL2=cspice_spkezr('392',ttd,'ECLIPJ2000','NONE','Sun');
    rr_SEL2=xx_SEL2(1:3,:)./LU;

    rrt=cspice_spkpos(targ,ttd,'ECLIPJ2000','NONE','Sun')./LU;

    L=length(u);
    
    figure
    subplot(3,3,[1,4,7])
    plot3(rr_SEL2(1,1),rr_SEL2(2,1),rr_SEL2(3,1),'ob')
    hold on
    plot3(rrt(1,end),rrt(2,end),rrt(3,end),'kx')
    plot3(0,0,0,'+k')
    plot3(rr_SEL2(1,:),rr_SEL2(2,:),rr_SEL2(3,:),'color',[.8 .8 .8],'LineWidth',0.1)
    plot3(rrt(1,:),rrt(2,:),rrt(3,:),'color',[.8 .8 .8],'LineWidth',0.1)
    for i=1:L-1
        line([zz(i,1), zz(i+1,1)],[zz(i,2), zz(i+1,2)],[zz(i,3), zz(i+1,3)],'color',[u(i) 0 1-u(i)],'linewidth',u(i)+1);
    end
    view([45, 15])    
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    grid on
    grid minor
    xlabel('$x [AU]$')
    ylabel('$y [AU]$')
    zlabel('$z [AU]$')
    legend('SEL2', 'AST', 'Sun','Location','northeast');

end