% Hodgkin-Huxley style model for a neuron of the Suprachiasmatic Nucleus
% Modified from DeWoskin et al, PNAS, 2015 to include a persistent sodium current
% Currently in press

function [t,v] = solve_dfp(R,S,experiment)

% Input parameters:
% R and S set time of day
% experiment sets the INaP parameters depending on the type of experiment
%   -1: application of GSK3 inhibitor CHIR
%    1: constitutively active GSK3 double knockin
%    2: application of INaP blocker riluzole

% State variables are stored in vector v:
% V = v(1)
% m = v(2)
% h = v(3)
% n = v(4)
% rl = v(5)
% rnl = v(6)
% fnl = v(7)
% p = v(8)
% s = v(9)
% cas = v(10)
% cac = v(11)


chir = 0;
dki = 0;
riluzole = 0;

if experiment == -1
    chir = 1;
elseif experiment == 1
    dki = 1;
elseif experiment == 2
    riluzole = 1;
end

% Timespan
tspan = [0 4000];

% Initial condition
v0 = [-68.019892548251022,0.017309397927220,0.037593630215027,0.328377697311672,0.004382806955149,0.001691138120461,0.043439785470156,0.061223626073877,0.078447944612010,0.000006199169984,0.000118905830122];

%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
	ks=1.65e-4; 
	ts=0.1;
	bs=0.0;
	kc=8.59e-9;
	tc=1.75e3;
	bc=3.1e-8;

	taurl=3.1;
	taurnl=3.1;

	K1=3.93e-5;
	K2=6.55e-4;

    c=5.7;
    gk=3;
    gcal=6;
    gcanl=20;
    gkca = 198.0/(1.0+exp(R))+2.0;
    gkleak = 0.2/(1.0+exp(R-1));
    gnaleak = 0.0576;
    gna=229;

    if dki
        gnap=2.13+0.14./(1.0+exp(-S));
    elseif chir
        gnap=1.46+.51./(1.0+exp(-S));
    elseif riluzole
        gnap=1;    
    else
        gnap = 1.59+0.5./(1.0+exp(-S));
    end
    
    Ena=45;
    Ek=-97;
    Eca=54;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the ODE solver ode15s.
[t,v] = ode15s(@df,tspan,v0);

% Calculation of individual currents:
% Ina=gna*v(:,2).^3.*v(:,3).*(v(:,1)-Ena);
% Ik=gk*v(:,4).^4.*(v(:,1)-Ek);
% Ical=gcal*v(:,5).*(K1./(K2+v(:,10))).*(v(:,1)-Eca);
% Icanl=gcanl*v(:,6).*v(:,7).*(v(:,1)-Eca);
% Ikca=gkca*v(:,9).^2.*(v(:,1)-Ek);
% Ikleak=gkleak*(v(:,1)-Ek);
% Inaleak=gnaleak*(v(:,1)-Ena);
% Inap=gnap*v(:,8).*(v(:,1)-Ena);

function dvdt = df(t,v)

    dvdt(1) = 1.0/c*(-gna*v(2)^3*v(3)*(v(1)-Ena)-gk*v(4)^4*(v(1)-Ek)-gcal*v(5)*(K1/(K2+v(10)))*(v(1)-Eca)-gcanl*v(6)*v(7)*(v(1)-Eca)-gkca*v(9)^2*(v(1)-Ek)-gkleak*(v(1)-Ek)-gnaleak*(v(1)-Ena)-gnap*v(8)*(v(1)-Ena));
    dvdt(2) = (minf(v(1))-v(2))/taum(v(1));
    dvdt(3) = (hinf(v(1))-v(3))/tauh(v(1));
    dvdt(4) = (ninf(v(1))-v(4))/taun(v(1));
    dvdt(5) = (rlinf(v(1))-v(5))/taurl;
    dvdt(6) = (rnlinf(v(1))-v(6))/taurnl;
    dvdt(7) = (fnlinf(v(1))-v(7))/taufnl(v(1));
    dvdt(8) = (pinf(v(1))-v(8))/taup(v(1));
    dvdt(9) = (sinf(v(10))-v(9))/taus(v(10));
    dvdt(10) = -ks*(gcal*v(5)*(K1/(K2+v(10)))*(v(1)-Eca)+gcanl*v(6)*v(7)*(v(1)-Eca))-v(10)/ts+bs;
    dvdt(11) = -kc*(gcal*v(5)*(K1/(K2+v(10)))*(v(1)-Eca)+gcanl*v(6)*v(7)*(v(1)-Eca))-v(11)/tc+bc;

    dvdt=dvdt';

    function y=minf(x)
        y=1/(1+exp(-(x+35.2)/8.1));
    end
    function y=hinf(x)
        y=1/(1+exp((x+62)/2));
    end
    function y=ninf(x)
        y=1/(1+exp((x-14)/(-17)))^.25;
    end
    function y=rlinf(x)
        y=1/(1+exp(-(x+36)/5.1));
    end
    function y=rnlinf(x)
        y=1/(1+exp(-(x+21.6)/6.7));
    end
    function y=fnlinf(x)
        y=1/(1+exp((x+260)/65));
    end
    function y=sinf(x)
        y=1e7*x^2/(1e7*x^2+5.6);
    end
    function y=pinf(x)
        y=1/(1+exp(-(x+25)/7.4))^1.5;
    end
    function y=taum(x)
        y=exp(-(x+286)/160);
    end
    function y=tauh(x)
        y=0.51+exp(-(x+26.6)/7.1);
    end
    function y=taun(x)
        y=exp(-(x-67)/68);
    end
    function y=taufnl(x)
        y=exp(-(x-444)/220);
    end
    function y=taus(x)
        y=500/(1e7*x^2+5.6);
    end
    function y=taup(x)
        y=100;
    end
    end
end







  