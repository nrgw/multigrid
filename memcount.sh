#!/bin/bash

sum=0

for i in $(seq 1 $1)
do
        data=$({ /bin/time -f "%M" $2 > /dev/null; } 2>&1)
        sum=$(echo $sum+$data | bc)
done

avg=$(echo "scale=2; $sum / $1" | bc)
echo "Memory Usage Average: $avg kibibytes"
