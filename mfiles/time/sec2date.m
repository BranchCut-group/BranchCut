function dateStr = sec2date(secs, refDateStr)
%
% Seconds from reference date to date string
%
% Creator John Merryman - sarmap
% Date: 22 Mar 2017
%
% Usage:
%   dateStr = sec2date(secs, [refDateStr])
%
% s        : (input) Seconds elapsed since reference date.
% refDate  : (input)  Reference date (default = '01-Jan-1970').
% dateStr  : (output)  Date string in SARscape format
%
% Examples:
%     s1 = date2sec('08-OCT-2015 08:32:42.063759');
%     d1 = sec2date(s1)
%        d1 = 08-Oct-2015 08:32:42.06375909
%
if nargin < 2
    refDateStr = '01-Jan-1990';
end

ddmmyyyy = datestr(floor(secs/86400) + datenum(refDateStr),'dd-mmm-yyyy');

sec = secs - 86400*(floor(secs/86400)); % seconds of day
% NOTE: there seems to be a rounding issue with datevec.
% decimalSeconds below gets rounded to 1e-05 s
%[~,~,~,h,m,decimalSeconds] = datevec(sec / 86400);
%hmsStr = sprintf('%02d:%02d:%010.7f', h, m, decimalSeconds);
hmsStr = sec2hms(sec);
dateStr = upper([ddmmyyyy ' ' hmsStr]);
