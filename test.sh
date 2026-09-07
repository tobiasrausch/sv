#!/bin/bash

## Run installation and pixi shell first

if [ -d data/lr/ ]
then
    cat README.md | awk '/^```bash$/,/^```$/  {print} {next}' | grep -v '^`' > run.sh
    LINE=`cat run.sh | grep -n "cd data/lr/" | cut -f 1 -d ":"`
    tail -n +${LINE} run.sh > run.sh.tmp
    mv run.sh.tmp run.sh
    chmod a+x run.sh
    ./run.sh
    if [ $? -eq 0 ]
    then
	echo "Test completed successfully!"
    fi
    rm run.sh
fi
