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
    input_bam_name=$(dx describe --name "$input_bam")
    input_bam_index_name=$(dx describe --name "$input_bam_index")
    svd_prefix_name=$(dx describe --name "$svd_prefix")


    # Create variable for output filename by removing .bam extension
    bam_prefix="${input_bam_name%.bam}"

    # Download sample, reference, and svd prefix files

    dx-download-all-inputs --parallel

    # Create output directories
    mkdir -p /home/dnanexus/out/selfSM
    mkdir -p /home/dnanexus/out/ancestry

    # Location of verifyBamID2 Dockerfile
    docker_file_id=

    #read the DNA Nexus api key as a variable
    API_KEY_wquotes=$(echo $DX_SECURITY_CONTEXT |  jq '.auth_token')
    API_KEY=$(echo "$API_KEY_wquotes" | sed 's/"//g')
    echo "$API_KEY"

    # Get Docker image from 001_Tools
    dx download $docker_file_id --auth "${API_KEY}"
    docker_file=$(dx describe ${docker_file_id} --name)
    DOCKERIMAGENAME=`tar xfO ${docker_file} manifest.json | sed -E 's/.*"RepoTags":\["?([^"]*)"?.*/\1/'`
    docker load < /home/dnanexus/"${docker_file}"


    # Run verifyBamID2 using Docker

    docker run -v /home/dnanexus:/home/dnanexus/ \
        --rm \
        $DOCKERIMAGENAME \
        VerifyBamID2 \
        --SVDPrefix /home/dnanexus/in/svd_prefix/"$svd_prefix_name" \
        --Reference /home/dnanexus/in/reference_fasta/"$reference_fasta_name" \
        --BamFile /home/dnanexus/in/input_bam/"$input_bam_name" \
        --Output /home/dnanexus/out/"$bam_prefix" \
        --Verbose

    # Move output files to output directories
    mv "$bam_prefix".selfSM /home/dnanexus/out/selfSM
    mv "$bam_prefix".Ancestry /home/dnanexus/out/ancestry

    # Upload output files to DNAnexus
    dx-upload-all-outputs

fi
