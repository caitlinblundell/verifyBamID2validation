# Generates commands to copy samples into 003 validation project and run verifyBAMID on them
# Usage: bash generate_dx_commands.sh <source_project>

# Get source project id from command line argument
SOURCE_PROJECT=$1
VALIDATION_PROJECT=project-J0yfQ400Jy1pJP08gFkGq6q7

# Count number of BAM files in source project
BAM_FILES=$(dx find data --path $SOURCE_PROJECT:/output/ --name "*.bam" --norecurse --brief)
echo "Number of BAM files found: $(echo "$BAM_FILES" | wc -l)"

# Move into validation project
dx select $VALIDATION_PROJECT

# Generate batch TSV
BATCH_FILE="verifyBAMID_batch.tsv"
dx generate_batch_inputs \
--path $SOURCE_PROJECT:/output/ \
-o $BATCH_FILE \
-i input_bam="*.bam" \
-i input_bam_index="*.bam.bai"

# Run verifyBAMID applet on each BAM file in validation project as batch run
dx run project-ByfFPz00jy1fk6PjpZ95F27J:applet-J5KPy2Q0Jy1YFBxGZJbxpjPY \
--name "verifyBAMID_batchrun" \
--destination $VALIDATION_PROJECT:/cp2_samples/ \
--batch-tsv $BATCH_FILE