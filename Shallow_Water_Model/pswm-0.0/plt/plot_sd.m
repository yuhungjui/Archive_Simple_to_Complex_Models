function plot_sd( runid, rundir )
% PLOT_SD  Plots the scalar diagnostics from a pswm run
%
% Usage:  plot_sd( runid, rundir )
%
% Input:
%     runid  run identifier [default:  see get_runid]
%     rundir directory where output files are [default:  ../run]
%
% Output:
%     Plots all scalar diagnostics from the file rundir/runid.out
%
% See also:  GETSD, PLTFLD, GETFLD

% Author:  Scott R. Fulton
%
% Revision history:
% 30 Dec 2005 original version (for cswm)
% 08 Aug 2007 updated for pswm

% get the scalar diagnostics
if (~exist('rundir','var') | isempty(rundir) ), rundir = '../run/'; end
if (~exist('runid' ,'var') | isempty(runid ) ), runid  = get_runid(rundir); end
[t,smin,smax,savg,fcor,beta,cval,lfld,labelf] = get_sd( runid, rundir );
nt = length(t);

% figure the time scale
if (t(nt)<=60)
    tlabel = 'time (seconds)';
elseif (t(nt)<=60*60)
    t = t/60;
    tlabel = 'time (minutes)';
elseif (t(nt)<=24*60*60)
    t = t/(60*60);
    tlabel = 'time (hours)';
else
    t = t/(24*60*60);
    tlabel = 'time (days)';
end

% draw plots of variables
clf
lmin = 'b:';
lmax = 'g--';
lavg = 'k-';

% divergence
dmin = smin(:,4)/fcor;
dmax = smax(:,4)/fcor; 
davg = savg(:,4)/fcor;
subplot(2,2,1);
plot(t,dmax,lmax,t,davg,lavg,t,dmin,lmin);
axis([0 t(nt) -Inf Inf]);
xlabel(tlabel);
ylabel('\delta/f');
title('divergence');
%legend('maximum','mean','minimum')

% relative vorticity
zmin = smin(:,5)/fcor;
zmax = smax(:,5)/fcor; 
zavg = savg(:,5)/fcor;
subplot(2,2,2);
plot(t,zmax,lmax,t,zavg,lavg,t,zmin,lmin);
axis([0 t(nt) -Inf Inf]);
xlabel(tlabel);
ylabel('\zeta/f');
title('relative vorticity');
%legend('maximum','mean','minimum')

% deviation height
hmin = smin(:,3)/cval;
hmax = smax(:,3)/cval;
havg = savg(:,3)/cval;
subplot(2,2,3);
plot(t,hmax,lmax,t,havg,lavg,t,hmin,lmin);
axis([0 t(nt) -Inf Inf]);
xlabel(tlabel);
ylabel('(h-H)/H');
title('deviation height');
%legend('maximum','mean','minimum')

% enstrophy
Zmin = smin(:,8)/fcor^2;
Zmax = smax(:,8)/fcor^2;
Zavg = savg(:,8)/fcor^2;
subplot(2,2,4);
plot(t,Zmax,lmax,t,Zavg,lavg,t,Zmin,lmin);
axis([0 t(nt) 0 Inf]);
xlabel(tlabel);
ylabel('Z/f^2');
title('enstrophy');
%legend('maximum','mean','minimum')

pause;

% draw plots of conserved quantities

clf;

% mass
h = 1+havg;
subplot(3,2,1);
plot(t,h,'k-');
axis([0 t(nt) 0 2]);
xlabel(tlabel);
ylabel('h/H');
title('mass');

subplot(3,2,2);
h0 = h(1);
plot(t,100*(h-h0)/h0,'k-');
axis([0 t(nt) -Inf Inf]);
xlabel(tlabel);
ylabel('(h-h_0)/h_0  (percent)');
title('mass');

% potential vorticity
q = savg(:,7)/fcor; 
subplot(3,2,3);
plot(t,q,'k-');
axis([0 t(nt) 0 2]);
xlabel(tlabel);
ylabel('q/f');
title('potential vorticity');

subplot(3,2,4);
q0 = q(1);
plot(t,100*(q-q0)/q0,'k-');
axis([0 t(nt) -Inf Inf]);
xlabel(tlabel);
ylabel('(q-q_0)/q_0  (percent)');
title('potential vorticity');

% energies
K = savg(:, 9)/cval^3; K0 = K(1); if K0<=0, K0 = 1; end; labK = 'b--';
A = savg(:,10)/cval^3; A0 = A(1); if A0<=0, A0 = 1; end; labA = 'g:';
E = K+A;               E0 = E(1); if E0<=0, E0 = 1; end; labE = 'k-';
subplot(3,2,5);
plot(t,K,labK,t,A,labA,t,E,labE);
axis([0 t(nt) 0 Inf]);
xlabel(tlabel);
ylabel('energy/c^3');
title('energies');
legend('kinetic','potential','total')

subplot(3,2,6);
plot(t,100*(K-K(1))/K0,labK,t,100*(A-A(1))/A0,labA,t,100*(E-E(1))/E0,labE);
axis([0 t(nt) -Inf Inf]);
xlabel(tlabel);
ylabel('(E-E_0)/E_0  (percent)');
title('energies');
legend('kinetic','potential','total')
