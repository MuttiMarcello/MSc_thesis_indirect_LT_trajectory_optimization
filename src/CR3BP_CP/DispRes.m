function [prob] = DispRes(prob,output)
% Marcello Mutti - marcello1.mutti@mail.polimi.it - Politecnico Di Milano - 9/3/25
% CR3BP constant parameters

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

    % dimensional units
    LU=3.84405000e5;    % [km]
    TU=3.75676967e5;    % [s]
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

    rvtf=prob.xf;

    % terminal target error
    df=zz(end,1:6).'-rvtf;

    % throttle switching function profile
    ep=prob.epsilon;
    Se=SwFun(tt,zz,prob);

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

        plot3D(zz,u,prob);
    
        fprintf('ToF: %.4f days, Final Mass: %.4f kg, Propellant Mass: %.4f kg\n\n',ToF,mf,mp)
    
        fprintf('Position error %.4e km\nVelocity error %.4e km/s\n',norm(df(1:3))*LU, norm(df(4:6))*LU/TU)

        fprintf('Max absolute Hamiltonian variation %.4e%%\n\n',abs((max(H)-min(H)))*100)
    

        %-throttle---------------------------------------------------------
        subplot(3,2,2)
        plot(ttd,u)
        axis tight
        ylim([0 1.1])
        ylabel('$u$')
        grid on
        grid minor

        %-switching-function-----------------------------------------------
        subplot(3,2,4)
        plot(ttd,Se)
        axis tight
        ylabel('$S$')
        grid on
        grid minor
            
        %-mass-------------------------------------------------------------
        subplot(3,2,6)
        plot(ttd,zz(:,7)*m0)
        axis tight
        ylabel('$m\,[kg]$')
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

function plot3D(zz,u,prob)

    L=length(u);
    
    figure
    subplot(3,2,[1,3,5])
    plot(zz(1,1),zz(1,2),'ob')
    hold on
    plot(zz(end,1),zz(end,2),'kx')
    plot(-prob.param(3),0,'k+')
    for i=1:L-1
        line([zz(i,1), zz(i+1,1)],[zz(i,2), zz(i+1,2)],'color',[u(i) 0 1-u(i)],'linewidth',u(i)+1);
    end  
    xlim(1.1*[min(zz(:,1)) max(zz(:,1))])
    ylim(1.1*[-max(abs(zz(:,2))) max(abs(zz(:,2)))])
    grid on
    grid minor
    xlabel('$x [LU]$')
    ylabel('$y [LU]$')
    legend('GTO', 'Halo','Earth','Location','northeast');

end