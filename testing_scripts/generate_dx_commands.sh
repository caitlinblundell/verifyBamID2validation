# Generates commands to copy samples into 003 validation project and run verifyBAMID on them
# Usage: bash generate_dx_commands.sh <source_project>

# get source project id from command line argument
SOURCE_PROJECT=$1
VALIDATION_PROJECT=project-J0yfQ400Jy1pJP08gFkGq6q7

# get list of BAM files in source project using find data
BAM_FILES=$(dx find data --path $SOURCE_PROJECT:/output --name "*.bam" --norecurse --brief)
echo "Number of BAM files found: $(echo "$BAM_FILES" | wc -l)"

# move to validation project
dx select $VALIDATION_PROJECT
dx mkdir /cp2_samples

# copy each BAM file to validation project with its corresponding index file, and add to batch file
for BAM_FILE in $BAM_FILES; do
    # Get filename from file id
    BAM_NAME=$(dx describe "$BAM_FILE" --name)

    # Copy files
    dx cp $SOURCE_PROJECT:$BAM_FILE $VALIDATION_PROJECT:/cp2_samples/
    dx cp $SOURCE_PROJECT:/output/${BAM_NAME}.bai $VALIDATION_PROJECT:/cp2_samples/

    # Create batch TSV file to run as batch job
    BATCH_FILE="verifyBAMID_batch.tsv"
    echo -e "input_bam\tinput_bam_index\tskip" > $BATCH_FILE

    # run verifyBAMID applet on each BAM file in validation project
    dx run project-ByfFPz00jy1fk6PjpZ95F27J:applet-J5KPy2Q0Jy1YFBxGZJbxpjPY --name "verifyBAMID_${BAM_NAME%.bam}" --destination $VALIDATION_PROJECT:/cp2_samples/ -iinput_bam=$BAM_NAME -iinput_bam_index=${BAM_NAME}.bai -iskip=false
done