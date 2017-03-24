awk '{if($2 == "Period") per = $4; if(NR < 6) print $0; else print $1*per, $2, $3, $4,$5,$6,$7}' $1.txt >$1_per.txt
