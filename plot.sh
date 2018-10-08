#!/bin/bash

gnuplot -persist > /dev/null 2>&1 << EOF
set logscale y
plot 'log.check' using 1:2 , 'log.check' using 1:3 

while (1){

replot
pause 1
}
EOF
