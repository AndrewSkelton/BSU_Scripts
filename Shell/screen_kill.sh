#!/bin/bash

for session in $(screen -ls | grep -o '[0-9]\{5\}')
do
	screen -S "${session}" -X quit;
done

for session in $(screen -ls | grep -o '[0-9]\{4\}')
do      
        screen -S "${session}" -X quit;
done

