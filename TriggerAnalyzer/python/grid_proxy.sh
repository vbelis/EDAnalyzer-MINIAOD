#!bin/bash
   echo -n "Want default time validity of the proxy (12h)? (y/n) "
   read   answer
   if [ $answer = "y"  ] 
      then
         voms-proxy-init --rfc --voms cms
   else 
       echo -n "How much validity time do you want in XY:ZK (h:m)?"
       read validity_time
   voms-proxy-init --rfc --valid $validity_time --voms cms

   fi 
