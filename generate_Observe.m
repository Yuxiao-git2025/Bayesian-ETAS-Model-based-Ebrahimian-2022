function [timein,magin,lonin,latin]=generate_Observe()
global T_start T_end time M Mc longitude latitude
ind=(time>=T_start & time<T_end & M>=Mc);
magin=M(ind);
magin=magin(:);

timein=time(ind);
timein=timein(:);

lonin=longitude(ind);
lonin=lonin(:);

latin=latitude(ind);
latin=latin(:);

end
