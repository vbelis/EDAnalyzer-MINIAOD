#!bin/bash
   echo -n "Do you need more that 1h of time validty for the proxy? (y/n) "
   read   answer
   if [ $answer = "n"  ] 
      then
         voms-proxy-init --rfc --valid 01:00 --voms cms
   else 
       echo -n "How much validity time do you want in XY:ZK (h:m)?"
       read validity_time
   voms-proxy-init --rfc --valid $validity_time --voms cms

   fi 
