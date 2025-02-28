% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% Procedure 4. Time Integration and Calculation, Output.
% 
% Integrate through time from Start Time to End Time and calculate the
% precipitation with convective mass flux.
% 
% =========================================================================
% 
% Time Integration:
% 
% =========================================================================


for i = 1:smax;
    
    Fso         = rho0.*CD.*abs(Vs).*(S00-SM(i));
    Fqo         = rho0.*CD.*abs(Vs).*(q00-qM(i));
    a           = (SM(i)+Lc.*qM(i)-SBt(i))./(PT-PB);
    
    % LHS         = (R-1).*(g.*(SBt(i)-SM(i))./del_PM+g.*a+g.*(1-c)./(PB-PT).*Lc.*qM(i))+(Lc.*R-W).*(g./del_PM.*(qM(i)-qBt(i)));
    LHS         = (1-R).*(g.*(SBt(i)-SM(i))./del_PM+g.*a+g.*(1-c).*Lc.*qM(i)./(PB-PT))+(W-Lc.*R).*(g./del_PM).*(qM(i)-qBt(i));
    % RHS         = (1-R).*(g.*Fso./del_PM-(SM(i)-SRM)./tao+(SBt(i)-SRB)./tao)-(W-Lc.*R).*(g.*Fqo./del_PM);
    RHS         = (R-1).*(g.*Fso./del_PM-(SM(i)-SRM)./tao+(SBt(i)-SRB)./tao)+(W-Lc.*R).*(g.*Fqo./del_PM);
    
    Mc(i)       = RHS./LHS; % Convective Mass Flux (kg/(m^2*s))
    
    if ( Mc(i) <= 0 )
        Mc(i) = 0;
    end
        
    dSM         = (g.*Fso./del_PM-(SM(i)-SRM)./tao+g.*Mc(i).*(SBt(i)-SM(i))./del_PM).*dt;
    dSBt        = (-(SBt(i)-SRB)./tao-g.*Mc(i).*a-g.*(1-c)./(PB-PT).*Lc.*qM(i).*Mc(i)).*dt;
    dqM         = (g.*Fqo/del_PM-g.*Mc(i)./del_PM.*(qM(i)-qBt(i))).*dt;
    dqBt        = ((1-c).*qM(i)-qBt(i)).*g.*Mc(i)./Pq.*dt;
    
    SM(i+1)     = SM(i) + dSM;
    SBt(i+1)    = SBt(i) + dSBt; 
    qM(i+1)     = qM(i) + dqM;
    qBt(i+1)    = qM(i) + dqBt;

end

clear i

% =========================================================================
% 
% Calculation:
% 
% =========================================================================

PR  = c.*qM(1:smax).*Mc(1:smax); 
    % precipitation assumed according to convective mass flux (kg/(m^2*s))
    
PR(PR<=0) = 0;

% =========================================================================
