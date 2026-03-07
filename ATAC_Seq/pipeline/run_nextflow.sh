#!/bin/bash
        
nextflow run nf-core/atacseq \
	-profile docker \
	-params-file ./params.yaml \
	-r 2.1.2 \
	-c config.config \
	-resume
