
#!/bin/bash

echo -n "File Name you want to xrdcp: "
read  file_name

xrdcp root://cms-xrd-global.cern.ch/$file_name . 
