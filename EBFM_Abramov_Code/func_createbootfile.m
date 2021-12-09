%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save output file for rebooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_createbootfile(A,io)
boot = struct;
boot.subZ = A.subZ;
boot.subW = A.subW;
boot.subD = A.subD;
boot.subS = A.subS;
boot.subT = A.subT;
boot.subSS = A.subSS;
boot.subSOIL = A.subSOIL;
boot.subTmean = A.subTmean;
boot.mbal = A.mbal;
boot.mbal_stake = A.mbal_stake;
boot.mbal_snow = A.mbal_snow;
boot.snowmass = A.snowmass;
boot.Tsurf = A.Tsurf;
boot.ys = A.ys;
boot.timelastsnow = A.timelastsnow;
boot.alb_snow = A.alb_snow;

cd(io.rebootdir);
if (io.writebootfile)
    save(io.bootfileout,'boot');
end
cd(io.homedir);

end







