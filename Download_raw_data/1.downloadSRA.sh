#!/bin/env bash
# Download SRA sequences 

# Set the path to the destination folder
CSV_FILE="path/to/samplesheet_.csv"
OUTPUT_DIR="path/to/raw_data"

# Create the destination folder if it doesn't exist
mkdir -p "$OUTPUT_DIR/LSFscripts"
mkdir -p "$OUTPUT_DIR/error"

cd "$OUTPUT_DIR" || exit

tail -n +2 "$CSV_FILE" | while IFS=',' read -r PATIENT SAMPLE LANE FASTQ_1 FASTQ_2; do
    SAMPLE_NAME="$SAMPLE"
    BATCH_SCRIPT="$OUTPUT_DIR/LSFscripts/${SAMPLE_NAME}_batch.sh"
    ERROR_LOG="$OUTPUT_DIR/error/${SAMPLE_NAME}.log"

    sleep 10

    cat <<EOT > "$BATCH_SCRIPT"
#!/bin/bash
    # Start the timer
    start_time=\$(date +%s.%N)

    # Download the files using fasterq-dump with 16 threads and Measure CPU time
    { time prefetch "${SAMPLE_NAME}" && fasterq-dump -e 16 -f -x -p "${SAMPLE_NAME}" -O "$OUTPUT_DIR" >/dev/null; } 2>&1

    # Capture the exit status of the previous command
    exit_status=\$?

    # Calculate the runtime
    end_time=\$(date +%s.%N)
    runtime=\$(echo "\$end_time - \$start_time" | bc)

    # Check if the SRA has enough 2 reads and are renamed files 
    if [[ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_1.fastq.gz" && -f "${OUTPUT_DIR}/${SAMPLE_NAME}_2.fastq.gz" ]]; then
    
        # Change the name of the downloaded files
        mv "${OUTPUT_DIR}/${SAMPLE_NAME}_1.fastq" "${OUTPUT_DIR}/${SAMPLE_NAME}_1.fastq"
        mv "${OUTPUT_DIR}/${SAMPLE_NAME}_2.fastq" "${OUTPUT_DIR}/${SAMPLE_NAME}_2.fastq"

        # Compress the downloaded files
        gzip "${OUTPUT_DIR}/${SAMPLE_NAME}_1.fastq"
        gzip "${OUTPUT_DIR}/${SAMPLE_NAME}_2.fastq"
        echo "Compression complete"

        # Print the metrics for the current files
        echo "Files: ${SAMPLE_NAME}_1.fastq.gz and ${SAMPLE_NAME}_2.fastq.gz"
        echo "Runtime: \$runtime seconds"

        # Exit with the captured exit status
        exit \$exit_status
    
    else
        echo "Error: Downloaded files not successful for ${SAMPLE_NAME}"
    fi
EOT

    chmod +x "$BATCH_SCRIPT"
    # Execute the batch script and redirect both stdout and stderr to the error log
    "$BATCH_SCRIPT" > "$ERROR_LOG" 2>&1

    # Capture the exit status of the batch script
    exit_status=$?

    # Check if the batch script encountered an error
    if [[ $exit_status -ne 0 ]]; then
        echo "Error: Batch script encountered an error for ${SAMPLE_NAME}"
    fi
done

