#!/bin/bash

# Set the path to the destination folder
destination_folder="$HOME/masterarbeit/raw_data"

# Create the destination folder if it doesn't exist
mkdir -p "$destination_folder/LSFscripts"
mkdir -p "$destination_folder/output_error"

cd "$destination_folder" || exit

while IFS=$'\t' read -r sra name_of_sample; do
    BATCHINFILE="$destination_folder/LSFscripts/${name_of_sample}_batch.sh"
    OUTPUT_ERROR="$destination_folder/output_error/${name_of_sample}_output_error.txt"

    sleep 10

    cat <<EOT > "$BATCHINFILE"
#!/bin/bash
    # Start the timer
    start_time=\$(date +%s.%N)

    # Download the files using fasterq-dump with 16 threads and Measure CPU time
    { time prefetch "$sra" && fasterq-dump -e 16 -f -x -p "$sra" -O "$destination_folder" >/dev/null; } 2>&1

    # Capture the exit status of the previous command
    exit_status=\$?

    # Calculate the runtime
    end_time=\$(date +%s.%N)
    runtime=\$(echo "\$end_time - \$start_time" | bc)

    #  Check if the SRA have enough 2 reads and are renamed files 
    if [[ -f "${destination_folder}/${sra}_1.fastq.gz" && -f "${destination_folder}/${sra}_2.fastq.gz" ]]; then
    
        # Change the name of the downloaded files
        mv "${destination_folder}/${sra}_1.fastq" "${destination_folder}/${name_of_sample}_1.fastq"
        mv "${destination_folder}/${sra}_2.fastq" "${destination_folder}/${name_of_sample}_2.fastq"

        # Compress the downloaded files
        gzip "${destination_folder}/${name_of_sample}_1.fastq"
        gzip "${destination_folder}/${name_of_sample}_2.fastq"
        echo "Compression complete"

        # Print the metrics for the current files
        echo "Files: ${name_of_sample}_1.fastq.gz and ${name_of_sample}_2.fastq.gz"
        echo "Runtime: \$runtime seconds"

        # Move the compressed files to the specified destination folder
        mv "${destination_folder}/${name_of_sample}_1.fastq.gz" "/mnt/sda/Masterarbeit/raw_data/${name_of_sample}_1.fastq.gz"
        mv "${destination_folder}/${name_of_sample}_2.fastq.gz" "/mnt/sda/Masterarbeit/raw_data/${name_of_sample}_2.fastq.gz"

        # Exit with the captured exit status
        exit \$exit_status
    
    else
    
        echo "Error: Downloaded files not successful for ${name_of_sample}"
    fi
EOT

    chmod +x "$BATCHINFILE"
    # Execute the batch script
    "$BATCHINFILE" > "$OUTPUT_ERROR" 2>&1

    # Capture the exit status of the batch script
    exit_status=\$?

    # Check if the batch script encountered an error
    if [[ \$exit_status -ne 0 ]]; then
        echo "Error: Batch script encountered an error for ${name_of_sample}"
    fi
  
done < "$HOME/masterarbeit/list_SRA_raw_data_new.tsv"

