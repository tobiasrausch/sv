#!/bin/bash

if [ -d mamba/ ]
then
    if [ -d data/lr/ ]
    then
	cat README.md | awk '/^```bash$/,/^```$/  {print} {next}' | grep -v '^`' | grep -v "^igv" > run.sh
	chmod a+x run.sh
	./run.sh
	if [ $? -eq 0 ]
	then
	    echo "Test completed successfully!"
	fi
	rm run.sh
    fi
fi
