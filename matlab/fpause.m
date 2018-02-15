function [ ] = fpause(message)
%FPAUSE Equivalent to fortran PAUSE
%  The PAUSE statement suspends execution, and waits for you to type: go
%  Any other input will terminate the program
  prompt = 'Do you want to continue? Y/N: ';
  resp = input(prompt,'s');
  if isempty(resp) || ~strcmpi(resp,'Y')
      error('roex:fpause',['User terminated: ' message]);
  end
  
end

