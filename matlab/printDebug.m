function printDebug( formatStr, varargin )
%PRINTPROGRESS Prints variables passed with default format
% To use printDebug Crete a global variable DEBUG = true in any workspace
% global DEBUG; DEBUG = true
% printDebug('a=%1.20E   b=%1.20E   e=%1.20E   d=%1.20E   x=%1.20E\n',a,b,e,d,x);
   global DEBUG
   if DEBUG
      fprintf(formatStr,varargin{:});
   end
end

