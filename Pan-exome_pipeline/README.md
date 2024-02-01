# Step 6

## Up/Downstream pangenome WDL

### Purpose

**Up/Downstream Pangenome WDL** Workflow encompasses four primary stages: 1) Pan-exome short read mapping, 2) Germline variant calling, 3) Annotation, and 4) Comprehensive reporting. This versatile workflow is designed for use in both pangenome and pan-exome analyses.

The provided command initializes the execution of this workflow using cromwell. The workflow is defined by the WDL script at `path/to/up-downstream_pangenome.wdl` and collaborates with `path/to/biotask_utils.wdl` to articulate detailed command-line instructions for each task.

For more detailed information, please refer to the `parameter_meta` section in the [Pan-exome_pipeline/up-downstream_pangenome.wdl](https://github.com/LuciaNhuNguyen/Masterarbeit/blob/main/Pan-exome_pipeline/up-downstream_pangenome.wdl).

### Prerequisites

Before executing the command, ensure the following prerequisites are met:
- `cromwell`, specified by the `JAR` file `cromwell-86.jar` is installed.
- The required input `JSON` file `SAMPLE_NAME_inputs.json` is prepared.
- An options JSON file `options.json` is available.
- The necessary data files, such as raw FASTQ files, pangenome files (`Pan-exome.giraffe.gbz`, `Pan-exome.dist`, `Pan-exome.min`), reference files, liftover scripts, and configuration files, are correctly specified in the `JSON` input.

For supplementary files supporting the WDL, please locate them in the [Pan-exome_pipeline](https://github.com/LuciaNhuNguyen/Masterarbeit/tree/main/Pan-exome_pipeline).  

### Usage

Execute the following command:

```bash
/usr/bin/time -v java -jar cromwell-86.jar run "path/to/up-downstream_pangenome.wdl" -i "SAMPLE_NAME_inputs.json" -o "options.json" 2>&1 | tee "SAMPLE_NAME.log"
```

![img](https://github.com/LuciaNhuNguyen/Masterarbeit/blob/main/Figures/Pan_pipeline.png)
Figure 3: Pipeline for WES germline small-variants analysis in the pan-exome.
