# This script extracts the FREEMIX value and sample name from the selfSM files produced in a verifyBamID run
# Use with generate_dx_commands.sh to validate the results of many verifyBamID runs
# Usage: bash retrieve_results.sh <validation_project>

# Retrieve project id
VALIDATION_PROJECT=$1
dx select $VALIDATION_PROJECT
# Define input and output paths
INPUT_FOLDER_PATH="/cp2_samples/QC/verifybamid/"
OUTPUT_DIR="/cp2_samples/results_summary/"
# Ensure output directory exists
dx mkdir -p $OUTPUT_DIR

# Take run id from first sample in output folder
FIRST_SAMPLE_NAME=$(dx ls $INPUT_FOLDER_PATH | head -n 1)
RUN_ID="${FIRST_SAMPLE_NAME%%_*}"

# Create output file
OUTPUT_FILE="${RUN_ID}_verifyBamID_results_summary.tsv"
echo -e "Sample_name\tFREEMIX" > "$OUTPUT_FILE"

# Get list of selfSM file ids in the output folder
dx find data --path $VALIDATION_PROJECT:$INPUT_FOLDER_PATH --name "*.selfSM" --brief | \
# Loop through each file id
while read -r file_id; do
    # then stream file content to stdout using cat
    dx cat "$file_id" | \
    # then awk to extract FREEMIX value and sample name:
    # On the second line (NR==2), print first column
    awk 'NR==2 {print $1 "\t" $7}' >> "$OUTPUT_FILE"
done

# Upload summary file to DNAnexus
dx upload "$OUTPUT_FILE" --path $OUTPUT_DIR