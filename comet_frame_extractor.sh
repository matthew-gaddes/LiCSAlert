echo "started"

# specify the target directory
target_directory="/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/volc-portal/processing_output/licsalert/"

# outdir
outdir="./"

# List all directories in the target directory
directories=$(find "$target_directory" -maxdepth 1 -type d)

# Iterate over each directory
for directory in $directories; do
    # Extract the directory name from the path
    dir_name=$(basename "$directory")
    
    # Run ls on each directory and redirect the output to a file
    ls "$directory" > "$outdir/${dir_name}_output.txt"
    
    echo "Listing for $dir_name is saved in ${dir_name}_output.txt"
done
