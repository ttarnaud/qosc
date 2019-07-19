function APf = APfr(APtimes,Tstart,Tend,PLOT)
APt = APtimes(APtimes>=Tstart&APtimes<=Tend); APdt = circshift(APt,-1,2)-APt;
APt = APt(1:end-1); APdt = APdt(1:end-1);
if PLOT
hold on; plot(APt,1./APdt); hold off;
end
APf = 1./APdt;
end