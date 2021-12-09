function [OUT,io] = func_writetofile(OUT,io,A,grid,t,time,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save model output to files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify variables to be written
if t==1
    io.varsout = [];
    io.unitsout = [];
    io.descout = [];
    io.typeout = [];
    if io.out_surface
        io.varsout{end+1} = 'gridzmask';        io.unitsout{end+1} = 'm a.s.l.';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Elevation';
        io.varsout{end+1} = 'gridz';            io.unitsout{end+1} = 'm a.s.l.';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Elevation'; 
        io.varsout{end+1} = 'smb';              io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Mass balance';
        io.varsout{end+1} = 'Tsurf';            io.unitsout{end+1} = 'K';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Surface temperature';
        io.varsout{end+1} = 'climT';            io.unitsout{end+1} = 'K';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Air temperature';
        io.varsout{end+1} = 'climP';            io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Precipitation';
        io.varsout{end+1} = 'climC';            io.unitsout{end+1} = 'fraction';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Cloud cover';
        io.varsout{end+1} = 'climRH';           io.unitsout{end+1} = 'fraction';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Relative humidity';
        io.varsout{end+1} = 'climWS';           io.unitsout{end+1} = 'm s-1';           io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Wind speed';
        io.varsout{end+1} = 'climPres';         io.unitsout{end+1} = 'Pa';              io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Air pressure';
        io.varsout{end+1} = 'climrain';         io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Rainfall';
        io.varsout{end+1} = 'climsnow';         io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Snowfall';
        io.varsout{end+1} = 'climTlapse';       io.unitsout{end+1} = 'K m^{-1}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Temperature lapse rate';
        io.varsout{end+1} = 'climPreslapse';    io.unitsout{end+1} = 'Pa m^{-1}';       io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Pressure exponential lapse rate';
        io.varsout{end+1} = 'climPotlapse';     io.unitsout{end+1} = 'K m^{-1}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Potential temperature lapse rate';
        io.varsout{end+1} = 'mbal_snow';        io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Snow mass balance';
        io.varsout{end+1} = 'mbal';             io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Cumulative mass balance';
        io.varsout{end+1} = 'mbal_stake';       io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Cumulative stake mass balance';        
        io.varsout{end+1} = 'melt';             io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Melt';
        io.varsout{end+1} = 'refr';             io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Refreezing';
        io.varsout{end+1} = 'refr_seasnl_snow'; io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Refreezing seasonal snow';
        io.varsout{end+1} = 'refr_intacc';      io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Internal accumulation';
        io.varsout{end+1} = 'runoff';           io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Runoff';
        io.varsout{end+1} = 'runoff_surf';      io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Surface runoff';
        io.varsout{end+1} = 'runoff_slush';     io.unitsout{end+1} = 'm w.e.';          io.typeout{end+1} = 'sum';      io.descout{end+1} = 'Slush runoff';
        io.varsout{end+1} = 'SWin';             io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Incoming SW radiation';
        io.varsout{end+1} = 'SWout';            io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Reflected SW radiation';
        io.varsout{end+1} = 'LWin';             io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Incoming LW radiation';
        io.varsout{end+1} = 'LWout';            io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Outgoing LW radiation';
        io.varsout{end+1} = 'SHF';              io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Sensible heat flux';
        io.varsout{end+1} = 'LHF';              io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Latent heat flux';
        io.varsout{end+1} = 'GHF';              io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Subsurface heat flux';
        io.varsout{end+1} = 'surfH';            io.unitsout{end+1} = 'm';               io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Surface height';
        io.varsout{end+1} = 'alb';              io.unitsout{end+1} = 'fraction';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Albedo';
        io.varsout{end+1} = 'climq';            io.unitsout{end+1} = 'kg kg^{-1}';      io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Specific humidity';    
        io.varsout{end+1} = 'climVP';           io.unitsout{end+1} = 'Pa';              io.typeout{end+1} = 'mean';     io.descout{end+1} = 'Vapour pressure';      
        io.varsout{end+1} = 'insol_shade';      io.unitsout{end+1} = '-';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'insol shade';
        io.varsout{end+1} = 'insol_TOA';        io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'insol TOA' ;
        io.varsout{end+1} = 'insol_TOAflat';    io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'insol TOAflat';
        io.varsout{end+1} = 'insol_TOAdir';     io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'insol TOAdir';
        io.varsout{end+1} = 'insol_TOAdiff';    io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'insol TOAdiff';   
        io.varsout{end+1} = 'insol_TOAshade';   io.unitsout{end+1} = 'W m^{-2}';        io.typeout{end+1} = 'mean';     io.descout{end+1} = 'insol TOAshade';
        io.varsout{end+1} = 't_cl';             io.unitsout{end+1} = '-';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'tau_cl';
        io.varsout{end+1} = 't_w';              io.unitsout{end+1} = '-';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'tau_w';
        io.varsout{end+1} = 't_a';              io.unitsout{end+1} = '-';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'tau_a';
        io.varsout{end+1} = 't_rg';             io.unitsout{end+1} = '-';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'tau_rg';
        io.varsout{end+1} = 't_eff';            io.unitsout{end+1} = '-';               io.typeout{end+1} = 'mean';     io.descout{end+1} = 'tau_eff';

    end
    if io.out_subsurface
        io.varsout{end+1} = 'subD';             io.unitsout{end+1} = 'kg m^{-3}';       io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Density';
        io.varsout{end+1} = 'subT';             io.unitsout{end+1} = 'K';               io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Temperature';         
        io.varsout{end+1} = 'subS';             io.unitsout{end+1} = 'mm w.e.';         io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Slush water content'; 
        io.varsout{end+1} = 'subW';             io.unitsout{end+1} = 'mm w.e.';         io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Irreducible water';   
        io.varsout{end+1} = 'subZ';             io.unitsout{end+1} = 'm';               io.typeout{end+1} = 'sample';   io.descout{end+1} = 'Layer thickness';
    end
end



%% Update OUT.TEMP with variables to be stored
fn = io.varsout;
for i=1:numel(fn)
    temp_long = eval(['A.' fn{i}]);
    if t==1
        OUT.TEMP.(fn{i}) = zeros(size(temp_long));
    end
    if strcmp(io.typeout{i},'sample')
        if mod(t+floor(io.freqout/2),io.freqout)==0
            OUT.TEMP.(fn{i}) = temp_long';
        end
    elseif strcmp(io.typeout{i},'mean')
        OUT.TEMP.(fn{i}) = OUT.TEMP.(fn{i}) + temp_long/io.freqout;
    elseif strcmp(io.typeout{i},'sum')
        OUT.TEMP.(fn{i}) = OUT.TEMP.(fn{i}) + temp_long;
    end
end

%% Save output to binary files
cd(io.outdir);

if t==1
    for i=1:length(fn)
        io.fid(i) = fopen(['OUT_' fn{i} '.bin'], 'w');
    end
end

if mod(t,io.freqout)==0
    cd(io.outdir);
    for i=1:length(fn)
        OUT.(fn{i}) = OUT.TEMP.(fn{i});
        fwrite(io.fid(i),OUT.(fn{i}),'real*4');
        OUT.TEMP.(fn{i}) = 0.0;
    end
end

if t==time.tn
    for i=1:length(fn)
        fclose(io.fid(i));
    end
end

cd(io.homedir);

if year(datetime(datevec(time.TCUR)))>year(datetime(datevec(time.TCUR-time.dt)))
    func_createbootfile(A,io);
end

%% Save run info to file
if t==time.tn
    runinfo.grid = grid;
    runinfo.time = time;
    runinfo.IOout = io;
    runinfo.Cout = C;
    cd(io.outdir);
    save(io.infofile,'-struct','runinfo');
    cd(io.homedir);
end

end