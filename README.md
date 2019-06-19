# CodFISH
**Cod**on **F**requency **I**ndexing **S**oftware for **H**IV

## Usage

### Convert BAM/SAM file to codfish/nucfish format file

1. Install Docker CE (https://docs.docker.com/install/).

2. Download script:
   For exporting "nucfish" format, download "sam2nucfish-docker":

   ```bash
   curl -sL https://raw.githubusercontent.com/hivdb/codfish/master/bin/sam2nucfish-docker -o sam2nucfish-docker
   chmod +x sam2nucfish-docker
   ```

   For exporting "codfish" format, download "sam2codfish-docker":

   ```bash
   curl -sL https://raw.githubusercontent.com/hivdb/codfish/master/bin/sam2codfish-docker -o sam2codfish-docker
   chmod +x sam2codfish-docker
   ```

3. Use following command to convert BAM/SAM files.

   nucfish:
   ```bash
   ./sam2nucfish-docker /path/to/folders/containing/bam/files
   ```

   codfish:
   ```bash
   ./sam2codfish-docker /path/to/folders/containing/bam/files
   ```

   The script will automatically find all files with extension of .bam or .sam and convert them.

#### Troubleshoots

1. Error "`Padding (BAM_CPAD, X) is currently not supported. Please implement. Sorry about that.`"

   If you are using Geneious to export the BAM file, uncheck the "export padded CIGARs" option.

2. Error "`fetch called on bamfile without index`"

   If you are using Geneious to export the BAM file, check the "export BAM index file" option with the file extension ".bai".
