`/usr/bin/time -v java -jar cromwell-86.jar run "path/to/up-downstream_pangenome.wdl" -i "${SAMPLE_NAME}_inputs.json" -o "options.json" 2>&1 | tee "SRR7890849.log"`
