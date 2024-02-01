`/usr/bin/time -v java -jar cromwell-86.jar run "path/to/up-downstream_pangenome.wdl" -i "SAMPLE_NAME_inputs.json" -o "options.json" 2>&1 | tee "SAMPLE_NAME.log"`
