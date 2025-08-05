#!/bin/bash
# verifyBamID2_v2.0.1 0.0.2

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


# If skip is not set to true, then proceed
if [ "$skip" != "true" ]; then
 
    # Log input variables
       
    echo "Value of reference_fasta: $reference_fasta"
    echo "Value of input_bam: $input_bam"
    echo "Value of input_bam_index: $input_bam_index"
    echo "Value of svd_prefix: $svd_prefix"

    # Get input filenames as strings

    reference_fasta_name=$(dx describe --name "$reference_fasta")
    reference_fasta_index_name=$(dx describe --name "$reference_fasta_index")
    input_bam_name=$(dx describe --name "$input_bam")
    input_bam_index_name=$(dx describe --name "$input_bam_index")

    # Set path to svd files
    svd_prefix_path="/home/dnanexus/svd/${svd_prefix}"

    # Create variable for output filename by removing .bam extension
    bam_prefix="${input_bam_name%.bam}"

    # Download sample, reference, and svd prefix files

    dx-download-all-inputs --parallel

    # Copy bai into the input_bam directory
    
    cp /home/dnanexus/in/input_bam_index/"$input_bam_index_name" /home/dnanexus/in/input_bam/"${input_bam_name}.bai"

    # Decompress reference fasta if it is gzipped
    if [[ "$reference_fasta_name" == *.gz ]]; then
        gunzip -c /home/dnanexus/in/reference_fasta/"$reference_fasta_name" > /home/dnanexus/in/reference_fasta/"${reference_fasta_name%.gz}"
        reference_fasta_name="${reference_fasta_name%.gz}"
    fi


    # Unpack fasta-index files
    mkdir -p /home/dnanexus/in/reference_fasta_index_extracted
    tar -xzvf /home/dnanexus/in/reference_fasta_index/"$reference_fasta_index_name" \
        -C /home/dnanexus/in/reference_fasta_index_extracted

    # Move index files into the same directory as the reference FASTA
    mv /home/dnanexus/in/reference_fasta_index_extracted/* /home/dnanexus/in/reference_fasta/
    
   ## Unpack the saved Docker image tarball into a .tar file
    gunzip -c /home/dnanexus/verifybamid2.tar.gz > /home/dnanexus/verifybamid2.tar

# Load the Docker image
    DOCKERIMAGENAME=$(
        docker load < /home/dnanexus/verifybamid2.tar \
        | grep "Loaded image:" \
        | cut -d' ' -f3
    )

    # Remove the .tar to save space
    rm /home/dnanexus/verifybamid2.tar

    # Create output directories
    mkdir -p /home/dnanexus/out/selfSM/output_selfSM/
    mkdir -p /home/dnanexus/out/ancestry/output_ancestry/
    mkdir -p /home/dnanexus/out/"$bam_prefix"

    # Run verifyBamID2

    docker run -v /home/dnanexus:/home/dnanexus/ \
        --rm \
        "$DOCKERIMAGENAME" \
        --SVDPrefix "$svd_prefix_path" \
        --Reference /home/dnanexus/in/reference_fasta/"$reference_fasta_name" \
        --BamFile /home/dnanexus/in/input_bam/"$input_bam_name" \
        --Output /home/dnanexus/out/"$bam_prefix" \
        --Verbose


    # Move output files to output directories
    mv /home/dnanexus/out/"$bam_prefix".selfSM /home/dnanexus/out/selfSM/output_selfSM/
    mv /home/dnanexus/out/"$bam_prefix".Ancestry /home/dnanexus/out/ancestry/output_ancestry/


    # Upload output files to DNAnexus
    dx-upload-all-outputs

fi
