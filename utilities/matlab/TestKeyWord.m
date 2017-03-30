function success = TestKeyWord(newLine,keyWord)

  if( length(keyWord)== 0)
    success = 0;
  else
    success = strncmp(newLine,keyWord,length(keyWord));
  end
  
  
  %if( ~success )
  % message=sprintf('%s%s%s','Sorry, key word=',keyWord,' not recognized. Exiting!');
  %  disp(message);
  %end
  %else
  %message=sprintf('%s%s%s','Key word=',keyWord,'recognized. Proceeding!');
  % disp(message);
  % return;
  %end;
  return;