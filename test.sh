#!/bin/bash

## Run installation, pixi shell and pixi download first

if [ -d data/lr/ ]
then
    cat README.md | awk '/^```bash$/,/^```$/  {print} {next}' | grep -v '^`' > run.sh
    LINE=`cat run.sh | grep -n "cd data/lr/" | cut -f 1 -d ":"`
    echo 'set -e' > run.sh.tmp
    tail -n +${LINE} run.sh >> run.sh.tmp
    mv run.sh.tmp run.sh
    chmod a+x run.sh
    ./run.sh
    rc=$?
    rm run.sh
    if [ ${rc} -eq 0 ]
    then
	echo "Test completed successfully!"
    else
	echo "Test FAILED (exit ${rc})"
	exit ${rc}
    fi
fi
