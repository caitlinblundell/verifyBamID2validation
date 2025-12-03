#!/bin/bash
# verifyBamID2_v2.0.1

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -o pipefail


# define log-step error handler that will log error to DNAnexus job log
log-step() {
    echo '{"error": {"type": "AppError", "message": "Error while '"$@"'; please refer to the job log for more details."}}' > ~/job_error.json
}

# If skip is not passed in, then proceed
if [ "$skip" != "true" ]; then


    log-step "downloading input variables"
    # Download sample, reference, and svd files
    dx-download-all-inputs --parallel

    # Log input variables     
    echo "Name of reference_fasta: $reference_fasta_name"
    echo "Name of input_bam: $input_bam_name"
    echo "Name of input_bam_index: $input_bam_index_name"

    # Copy bai into the input_bam directory to make it accessible to verifyBamID2
    cp /home/dnanexus/in/input_bam_index/"$input_bam_index_name" /home/dnanexus/in/input_bam/"${input_bam_name}.bai"


    log-step "validating svd files"
    # Get SVD prefix from first file in array
    svd_prefix="${svd_array_prefix[0]}" # helper bash variable _prefix removes suffixes defined in dxapp.json

    # Move all SVD files into a single directory as VBID2 expects
    mkdir -p /home/dnanexus/in/collated_SVD_files
    find /home/dnanexus/in/svd_array/ -type f -exec mv {} /home/dnanexus/in/collated_SVD_files/ \;

    # Validate all required SVD files exist
    missing_files=()
    for ext in bed mu UD V; do
        if [[ ! -f "/home/dnanexus/in/collated_SVD_files/${svd_prefix}.${ext}" ]]; then
        missing_files+=("${svd_prefix}.${ext}")
        fi
    done

    if [[ ${#missing_files[@]} -gt 0 ]]; then
        echo "ERROR: Missing SVD files: ${missing_files[*]}"
        exit 1
    fi


    log-step "validating reference genome"
    if  [[ "$reference_fasta_name" == *.tar* ]] # error if reference is a tarball
	    then
            echo "Received file: $reference_fasta_name"
		    echo "Invalid format: reference genome file is a tarball, provide .fa or .fa.gz"
		    exit 1
    fi
    # Decompress reference fasta if it is gzipped
    if [[ "$reference_fasta_name" == *.gz ]] 
        then 
            gunzip -c /home/dnanexus/in/reference_fasta/"$reference_fasta_name" > /home/dnanexus/in/reference_fasta/"${reference_fasta_name%.gz}"
            reference_fasta_name="${reference_fasta_name%.gz}"
            echo "Reference genome decompressed"
    fi

    # Check and extract reference fasta index:
    # Check reference index file tarball (.tar or .tar.gz)
    if [[ "$reference_fasta_index_name" != *.tar* ]]; then
        echo "Received file: $reference_fasta_index_name"
        echo "Invalid format: reference fasta .fai and .gzi should be packaged in a tarball"
        exit 1
    fi

    # If gzipped, decompress first
    if [[ "$reference_fasta_index_name" == *.gz ]]; then
        echo "Reference index is gzipped, decompressing..."
        gunzip -c "/home/dnanexus/in/reference_fasta_index/$reference_fasta_index_name" \
            > "/home/dnanexus/in/reference_fasta_index/${reference_fasta_index_name%.gz}"
        reference_fasta_index_name="${reference_fasta_index_name%.gz}"
        echo "Reference index decompressed to $reference_fasta_index_name"
    fi

    # Extract reference index tarball
    mkdir -p /home/dnanexus/in/reference_fasta_index_extracted
    tar -xvf "/home/dnanexus/in/reference_fasta_index/$reference_fasta_index_name" \
        -C /home/dnanexus/in/reference_fasta_index_extracted

    # Move index files into the same directory as the reference FASTA
    mv /home/dnanexus/in/reference_fasta_index_extracted/* /home/dnanexus/in/reference_fasta/


    log-step "preparing Docker image"
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

    # Create variable for output filename by removing .bam extension
    bam_prefix="${input_bam_name%.bam}"

    # Create output directories
    mkdir -p /home/dnanexus/out/all_outputs/vbid_output
    mkdir -p /home/dnanexus/out/"$bam_prefix"


    log-step "running verifyBamID2"
    # Temporarily disable exit-on-error so that files (e.g. pileup) can be examined if VBID2 fails
    set +e
    # Run verifyBamID2
    docker run -v /home/dnanexus:/home/dnanexus/ \
        --rm \
        "$DOCKERIMAGENAME" \
        --SVDPrefix /home/dnanexus/in/collated_SVD_files/"$svd_prefix" \
        --Reference /home/dnanexus/in/reference_fasta/"$reference_fasta_name" \
        --BamFile /home/dnanexus/in/input_bam/"$input_bam_name" \
        --Output /home/dnanexus/out/"$bam_prefix" \
        $additional_flags
    status=$?
    set -e  # re-enable exit-on-error

    log-step "uploading output files"

    missing_outputs=()

    # Handle .selfSM, .Ancestry, and .Pileup consistently
    for ext in selfSM Ancestry Pileup; do
        src1="/home/dnanexus/out/${bam_prefix}.${ext}"
        src2="/home/dnanexus/out/${bam_prefix}/${bam_prefix}.${ext}"
        if [[ -f "$src1" ]]; then
            mv "$src1" /home/dnanexus/out/all_outputs/vbid_output/
        elif [[ -f "$src2" ]]; then
            mv "$src2" /home/dnanexus/out/all_outputs/vbid_output/
        else
            if [[ "$ext" == "selfSM" || "$ext" == "Ancestry" ]]; then
                missing_outputs+=("${bam_prefix}.${ext}")
                echo "DEBUG: Missing $ext output"
            fi
        fi
    done

    if [[ ${#missing_outputs[@]} -gt 0 ]]; then
        echo "WARNING: One or both output files (.selfSM, .Ancestry) are missing: ${missing_outputs[*]}"
        echo "This can happen when too few SNP markers overlap between BAM and SVD panel."
    fi

    dx-upload-all-outputs

    if [[ $status -ne 0 ]]; then
        echo "WARNING: verifyBamID2 failed with exit code $status. Partial outputs (if any) have been uploaded."
    fi
    exit $status


fi
